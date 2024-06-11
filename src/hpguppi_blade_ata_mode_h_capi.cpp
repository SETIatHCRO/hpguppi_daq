#include <cassert>
#include <memory>

#include "blade/plan.hh"
#include "blade/base.hh"
#include "blade/logger.hh"
#include "blade/runner.hh"
#include "blade/pipelines/ata/mode_b.hh"
#include "blade/pipelines/generic/mode_h.hh"

extern "C" {
#include "hpguppi_blade_ata_mode_h_capi.h"
}

using namespace Blade;
using namespace Blade::Pipelines::ATA;
using namespace Blade::Pipelines::Generic;

using BladePipelineB = ModeB<CF32>;
using BladePipelineH = ModeH<CF32, F32>;

static Vector<Device::CPU, F64> blockJulianDate({1});
static Vector<Device::CPU, F64> blockDut1({1});

static struct {
    U64 StepCount = 0;
    void* UserData = nullptr;
    std::unordered_map<U64, void*> InputPointerMap;
    std::unordered_map<U64, void*> OutputPointerMap;
    std::unordered_map<U64, size_t> InputIdMap;
    std::unordered_map<U64, size_t> OutputIdMap;    

    struct {
        std::unique_ptr<Runner<BladePipelineB>> B; 
        std::unique_ptr<Runner<BladePipelineH>> H; 
    } RunnersInstances;

    struct {
        blade_stateful_cb* InputBufferPrefetch;
        blade_input_buffer_fetch_cb* InputBufferFetch;
        blade_input_buffer_enqueued_cb* InputBufferEnqueued;
        blade_input_buffer_ready_cb* InputBufferReady;
        blade_output_buffer_fetch_cb* OutputBufferFetch;
        blade_output_buffer_ready_cb* OutputBufferReady;

        blade_clear_queued_cb* InputClear;
        blade_clear_queued_cb* OutputClear;
    } Callbacks;
} State;


bool blade_ata_h_initialize(
    struct blade_ata_mode_h_config ata_h_config,
    size_t numberOfWorkers,
    struct blade_ata_observation_meta* observationMeta,
    struct LonLatAlt* arrayReferencePosition,
    double* obs_phase_center_radecrad,
    double* beamCoordinates_radecrad,
    double* antennaPositions_xyz,
    double _Complex* antennaCalibrations
) {
    if (State.RunnersInstances.B || State.RunnersInstances.H) {
        BL_FATAL("Can't initialize because Blade Runner is already initialized.");
        throw Result::ASSERTION_ERROR;
    }

    std::vector<XYZ> antennaPositions(ata_h_config.inputDims.NANTS);
    std::vector<RA_DEC> beamCoordinates(ata_h_config.beamformerBeams);
    std::vector<std::complex<double>> antennaCalibrationsCpp(
            ata_h_config.inputDims.NANTS*\
            ata_h_config.inputDims.NCHANS*\
            ata_h_config.inputDims.NPOLS);
    int i;
    for(i = 0; i < ata_h_config.inputDims.NANTS; i++){
        antennaPositions[i].X = antennaPositions_xyz[i*3 + 0];
        antennaPositions[i].Y = antennaPositions_xyz[i*3 + 1];
        antennaPositions[i].Z = antennaPositions_xyz[i*3 + 2];
    }
    for(i = 0; i < ata_h_config.beamformerBeams; i++){
        beamCoordinates[i].RA = beamCoordinates_radecrad[i*2 + 0];
        beamCoordinates[i].DEC = beamCoordinates_radecrad[i*2 + 1];
    }
    memcpy(antennaCalibrationsCpp.data(), antennaCalibrations,
            antennaCalibrationsCpp.size()*sizeof(antennaCalibrationsCpp[0]));

    ArrayDimensions beamformerInputDimensions = ArrayDimensions({
        .A = ata_h_config.inputDims.NANTS,
        .F = ata_h_config.inputDims.NCHANS,
        .T = ata_h_config.inputDims.NTIME,
        .P = ata_h_config.inputDims.NPOLS,
    });

    ArrayTensor phasorAntennaCalibrations = ArrayTensor<Device::CPU, CF64>({
        ata_h_config.inputDims.NANTS,
        ata_h_config.inputDims.NCHANS * ata_h_config.channelizerRate,
        1,
        ata_h_config.inputDims.NPOLS,
    });

    const size_t calAntStride = 1;
    const size_t calPolStride = beamformerInputDimensions.numberOfAspects() * calAntStride;
    const size_t calChnStride = beamformerInputDimensions.numberOfPolarizations() * calPolStride;

    const size_t weightsPolStride = 1;
    const size_t weightsChnStride = beamformerInputDimensions.numberOfPolarizations() * weightsPolStride;
    const size_t weightsAntStride = beamformerInputDimensions.numberOfFrequencyChannels() * ata_h_config.channelizerRate * weightsChnStride;
    BL_INFO("Expanding the {} coarse-channel coefficients by a factor of {}.", beamformerInputDimensions.numberOfFrequencyChannels(), ata_h_config.channelizerRate);

    U64 inputIdx, frqIdx, outputIdx, antIdx, chnIdx, polIdx, fchIdx;
    for (antIdx = 0; antIdx < beamformerInputDimensions.numberOfAspects(); antIdx++) {
        for (chnIdx = 0; chnIdx < beamformerInputDimensions.numberOfFrequencyChannels(); chnIdx++) {
            for (polIdx = 0; polIdx < beamformerInputDimensions.numberOfPolarizations(); polIdx++) {
                inputIdx = chnIdx * calChnStride +
                    polIdx * calPolStride + 
                    antIdx * calAntStride;
                for (fchIdx = 0; fchIdx < ata_h_config.channelizerRate; fchIdx++) {
                    frqIdx = chnIdx * ata_h_config.channelizerRate + fchIdx;
                    outputIdx = antIdx * weightsAntStride +
                        polIdx * weightsPolStride +
                        frqIdx * weightsChnStride;

                    phasorAntennaCalibrations[outputIdx] = antennaCalibrationsCpp[inputIdx];
                }
            }
        }
    }

    State.RunnersInstances.B = Runner<BladePipelineB>::New(
        numberOfWorkers,
        {  
            .inputDimensions = beamformerInputDimensions,

            .preBeamformerChannelizerRate = ata_h_config.channelizerRate,
            // .preBeamformerPolarizerConvertToCircular = BLADE_ATA_MODE_H_CIRCULAR_POLARIZATION,

            .phasorObservationFrequencyHz = observationMeta->rfFrequencyHz,
            .phasorChannelBandwidthHz = observationMeta->channelBandwidthHz,
            .phasorTotalBandwidthHz = observationMeta->totalBandwidthHz,
            .phasorFrequencyStartIndex = observationMeta->frequencyStartIndex,
            .phasorReferenceAntennaIndex = observationMeta->referenceAntennaIndex,
            .phasorArrayReferencePosition = {
                .LON = arrayReferencePosition->LON,
                .LAT = arrayReferencePosition->LAT,
                .ALT = arrayReferencePosition->ALT
            },
            .phasorBoresightCoordinate = {
                .RA = obs_phase_center_radecrad[0],
                .DEC = obs_phase_center_radecrad[1]
            },
            .phasorAntennaPositions = antennaPositions,
            .phasorAntennaCalibrations = phasorAntennaCalibrations,
            .phasorBeamCoordinates = beamCoordinates,

            .beamformerIncoherentBeam = false,

            .detectorEnable = false,
            .detectorIntegrationSize = 1,
            .detectorNumberOfOutputPolarizations = 1,

            .castBlockSize = ata_h_config.castBlockSize,
            .channelizerBlockSize = ata_h_config.channelizerBlockSize,
            .phasorBlockSize = 512,
            .beamformerBlockSize = ata_h_config.beamformerBlockSize,
            .detectorBlockSize = 512,
        }
    );

    State.RunnersInstances.H = Runner<BladePipelineH>::New(
        numberOfWorkers, 
        {
            .inputDimensions = State.RunnersInstances.B->getWorker().getOutputBuffer().dims(),

            .accumulateRate = BLADE_ATA_MODE_H_ACCUMULATE_RATE, // ata_h_config.accumulateRate

            .detectorIntegrationSize = BLADE_ATA_MODE_H_INTEGRATION_SIZE, // ata_h_config.integrationSize
            .detectorNumberOfOutputPolarizations = BLADE_ATA_MODE_H_OUTPUT_NPOL,
        }
    );

    State.InputPointerMap.reserve(numberOfWorkers);
    State.OutputPointerMap.reserve(numberOfWorkers);
    State.InputIdMap.reserve(numberOfWorkers);
    State.OutputIdMap.reserve(numberOfWorkers);

    return true;
}

void blade_ata_h_terminate() {
    if (!State.RunnersInstances.B || !State.RunnersInstances.H) {
        BL_WARN("Can't terminate because Blade Runner isn't initialized.");
        // throw Result::ASSERTION_ERROR;
        // return;
    }

    for(const auto &[key, value]: State.InputIdMap) {
        State.Callbacks.InputClear(State.UserData, value);
    }
    State.InputIdMap.clear();
    for(const auto &[key, value]: State.OutputIdMap) {
        State.Callbacks.OutputClear(State.UserData, value);
    }
    State.OutputIdMap.clear();

    State.RunnersInstances.B.reset();
    State.RunnersInstances.H.reset();
}

U64 blade_ata_h_get_input_size() {
    assert(State.RunnersInstances.B);
    return State.RunnersInstances.B->getWorker().getInputBuffer().dims().size();
}

U64 blade_ata_h_get_output_size() {
    assert(State.RunnersInstances.H);
    return State.RunnersInstances.H->getWorker().getOutputBuffer().dims().size();
}

void blade_ata_h_set_block_time_mjd(double mjd) {
    blockJulianDate.data()[0] = mjd;
}

void blade_ata_h_set_block_dut1(double dut1) {
    blockDut1.data()[0] = dut1;
}

void blade_ata_h_register_user_data(void* user_data) {
    State.UserData = user_data;
}

void blade_ata_h_register_input_buffer_prefetch_cb(blade_stateful_cb* f) {
    State.Callbacks.InputBufferPrefetch = f;
}

void blade_ata_h_register_input_buffer_fetch_cb(blade_input_buffer_fetch_cb* f) {
    State.Callbacks.InputBufferFetch = f;
}

void blade_ata_h_register_input_buffer_enqueued_cb(blade_input_buffer_enqueued_cb* f) {
    State.Callbacks.InputBufferEnqueued = f;
}

void blade_ata_h_register_input_buffer_ready_cb(blade_input_buffer_ready_cb* f) {
    State.Callbacks.InputBufferReady = f;
}

void blade_ata_h_register_output_buffer_fetch_cb(blade_output_buffer_fetch_cb* f) {
    State.Callbacks.OutputBufferFetch = f;
}

void blade_ata_h_register_output_buffer_ready_cb(blade_output_buffer_ready_cb* f) {
    State.Callbacks.OutputBufferReady = f;
}

void blade_ata_h_register_blade_queued_input_clear_cb(blade_clear_queued_cb* f) {
    State.Callbacks.InputClear = f;
}

void blade_ata_h_register_blade_queued_output_clear_cb(blade_clear_queued_cb* f) {
    State.Callbacks.OutputClear = f;
}

size_t blade_ata_h_accumulator_counter() {
    return State.RunnersInstances.H->getNextWorker().getCurrentAccumulatorStep();
}

bool blade_ata_h_compute_step() {
    bool prefetch = State.Callbacks.InputBufferPrefetch(State.UserData);
    
    if(!State.RunnersInstances.B) {
        return false;
    }
    if(!State.RunnersInstances.H) {
        return false;
    }

    U64 callbackStep = 0;
    void* externalBuffer = nullptr;

    auto& ModeB = State.RunnersInstances.B; 
    auto& ModeH = State.RunnersInstances.H; 
    if (prefetch) {
        ModeB->enqueue([&](auto& worker) {
            // Check if next runner has free slot.
            Plan::Available(ModeH);

            size_t bufferId;
            // Calls client callback to request empty input buffer.
            if (!State.Callbacks.InputBufferFetch(State.UserData, &externalBuffer, &bufferId)) {
                Plan::Skip();
            }

            // Keeps track of pointer for "ready" callback.
            State.InputPointerMap.insert({State.StepCount, externalBuffer});
            State.InputIdMap.insert({State.StepCount, bufferId});

            // Create Memory::ArrayTensor from RAW pointer.
            auto input = ArrayTensor<Device::CPU, CI8>(externalBuffer, worker.getInputBuffer().dims());

            // Transfer input memory to the pipeline.
            Plan::TransferIn(worker, 
                            blockJulianDate,
                            blockDut1,
                            input);

            // Compute input data.
            Plan::Compute(worker);

            // Concatenate output data inside next pipeline input buffer.
            Plan::Accumulate(ModeH, ModeB, worker.getOutputBuffer());

            // Asynchronous CPU work
            State.Callbacks.InputBufferEnqueued(State.UserData, bufferId, ~0); // not optimal, move to spin-loop

            // Return job identity and increment counter.
            return State.StepCount++; 
        });
    }

    ModeH->enqueue([&](auto& worker) {
        // Try dequeue job from last runner. If unlucky, return.
        Plan::Dequeue(ModeB, &callbackStep);

        // If dequeue successfull, recycle input buffer.
        const auto& recycleBuffer = State.InputPointerMap[callbackStep];
        const auto& recycleBufferId = State.InputIdMap[callbackStep];
        State.Callbacks.InputBufferReady(State.UserData, recycleBuffer, recycleBufferId);
        State.InputPointerMap.erase(callbackStep);
        State.InputIdMap.erase(callbackStep);

        // Compute input data.
        Plan::Compute(worker);
        
        size_t outputBufferId;
        // Calls client callback to request empty output buffer.
        if (!State.Callbacks.OutputBufferFetch(State.UserData, &externalBuffer, &outputBufferId)) {
            Plan::Skip();
        }

        // Keeps track of pointer for "ready" callback.
        State.OutputPointerMap.insert({callbackStep, externalBuffer});
        State.OutputIdMap.insert({callbackStep, outputBufferId});

        // Create Memory::ArrayTensor from RAW pointer.
        auto output = ArrayTensor<Device::CPU, F32>(externalBuffer, worker.getOutputBuffer().dims());

        // Copy worker output to external output buffer.
        Plan::TransferOut(output, worker.getOutputBuffer(), worker);
        // Plan::TransferOutAsATPFrev(output, worker.getOutputBuffer(), worker);

        // Return job identity.
        return callbackStep;
    });

    // Dequeue last runner job and recycle output buffer.
    Result result = ModeH->dequeue(&callbackStep);
    if (result == Result::SUCCESS) {
    // if (ModeH->dequeue(&callbackStep)) {
        const auto& recycleBuffer = State.OutputPointerMap[callbackStep];
        const auto& recycleBufferId = State.OutputIdMap[callbackStep];
        State.Callbacks.OutputBufferReady(State.UserData, recycleBuffer, recycleBufferId);
        State.OutputPointerMap.erase(callbackStep);
        State.OutputIdMap.erase(callbackStep);
    }
    else if(result == Result::EXHAUSTED) {}
    else {
        BL_WARN("Result {}, was empty {}", result, ModeH->empty());
        // assert(0);
    }

    // Prevent memory clobber inside spin-loop.
    Plan::Loop();

    // Return if pipeline is computing something.
    return !(ModeB->empty() && ModeH->empty());
}
