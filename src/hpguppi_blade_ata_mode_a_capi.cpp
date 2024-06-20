#include <cassert>
#include <memory>

#include "blade/base.hh"
#include "blade/runner.hh"
#include "blade/bundles/ata/mode_b.hh"

extern "C" {
#include "hpguppi_blade_ata_mode_a_capi.h"
}

using namespace Blade;
using namespace Blade::Bundles::ATA;

using BladePipeline = ModeB<CI8, BLADE_ATA_MODE_A_OUTPUT_ELEMENT_T>;

static Tensor<Device::CPU, F64> blockJulianDate({1});
static Tensor<Device::CPU, F64> blockDut1({1});

static struct {
    U64 StepCount = 0;
    void* UserData = nullptr;
    std::unordered_map<U64, void*> InputPointerMap;
    std::unordered_map<U64, void*> OutputPointerMap;
    
    std::shared_ptr<Runner> pipelineRunner;
    std::shared_ptr<BladePipeline> pipeline;
    Duet<ArrayTensor<Device::CUDA, CI8>> inputBuffer;
    ArrayShape inputShape;
    ArrayShape outputShape;

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

bool blade_ata_a_initialize(
    struct blade_ata_mode_a_config ata_a_config,
    size_t numberOfWorkers,
    struct blade_ata_observation_meta* observationMeta,
    struct LonLatAlt* arrayReferencePosition,
    double* obs_phase_center_radecrad,
    double* beamCoordinates_radecrad,
    double* antennaPositions_xyz,
    double _Complex* antennaCalibrations
) {
    if (State.pipelineRunner) {
        BL_FATAL("Can't initialize because Blade Runner is already initialized.");
        throw Result::ASSERTION_ERROR;
    }

    std::vector<XYZ> antennaPositions(ata_a_config.inputDims.NANTS);
    std::vector<RA_DEC> beamCoordinates(ata_a_config.beamformerBeams);
    std::vector<std::complex<double>> antennaCalibrationsCpp(
            ata_a_config.inputDims.NANTS*\
            ata_a_config.inputDims.NCHANS*\
            ata_a_config.inputDims.NPOLS);
    int i;
    for(i = 0; i < ata_a_config.inputDims.NANTS; i++){
        antennaPositions[i].X = antennaPositions_xyz[i*3 + 0];
        antennaPositions[i].Y = antennaPositions_xyz[i*3 + 1];
        antennaPositions[i].Z = antennaPositions_xyz[i*3 + 2];
    }
    for(i = 0; i < ata_a_config.beamformerBeams; i++){
        beamCoordinates[i].RA = beamCoordinates_radecrad[i*2 + 0];
        beamCoordinates[i].DEC = beamCoordinates_radecrad[i*2 + 1];
    }
    memcpy(antennaCalibrationsCpp.data(), antennaCalibrations,
            antennaCalibrationsCpp.size()*sizeof(antennaCalibrationsCpp[0]));

    State.inputShape = ArrayShape({
        ata_a_config.inputDims.NANTS,
        ata_a_config.inputDims.NCHANS,
        ata_a_config.inputDims.NTIME,
        ata_a_config.inputDims.NPOLS,
    });
    
    State.outputShape = ArrayShape({
        ata_a_config.beamformerBeams + (BLADE_ATA_MODE_A_OUTPUT_INCOHERENT_BEAM ? 1 : 0),
        ata_a_config.inputDims.NCHANS * ata_a_config.channelizerRate,
        ata_a_config.inputDims.NTIME / (ata_a_config.channelizerRate * ata_a_config.integrationSize),
        ata_a_config.numberOfOutputPolarizations, // detector enabled
    });

    auto phasorAntennaCalibrations = ArrayTensor<Device::CPU, CF64>({
        ata_a_config.inputDims.NANTS,
        ata_a_config.inputDims.NCHANS * ata_a_config.channelizerRate,
        1,
        ata_a_config.inputDims.NPOLS,
    });

    const size_t calAntStride = 1;
    const size_t calPolStride = State.inputShape.numberOfAspects() * calAntStride;
    const size_t calChnStride = State.inputShape.numberOfPolarizations() * calPolStride;

    const size_t weightsPolStride = 1;
    const size_t weightsChnStride = State.inputShape.numberOfPolarizations() * weightsPolStride;
    const size_t weightsAntStride = State.inputShape.numberOfFrequencyChannels() * ata_a_config.channelizerRate * weightsChnStride;
    BL_INFO("Expanding the {} coarse-channel coefficients by a factor of {}.", State.inputShape.numberOfFrequencyChannels(), ata_a_config.channelizerRate);

    U64 inputIdx, frqIdx, outputIdx, antIdx, chnIdx, polIdx, fchIdx;
    for (antIdx = 0; antIdx < State.inputShape.numberOfAspects(); antIdx++) {
        for (chnIdx = 0; chnIdx < State.inputShape.numberOfFrequencyChannels(); chnIdx++) {
            for (polIdx = 0; polIdx < State.inputShape.numberOfPolarizations(); polIdx++) {
                inputIdx = chnIdx * calChnStride +
                    polIdx * calPolStride + 
                    antIdx * calAntStride;
                for (fchIdx = 0; fchIdx < ata_a_config.channelizerRate; fchIdx++) {
                    frqIdx = chnIdx * ata_a_config.channelizerRate + fchIdx;
                    outputIdx = antIdx * weightsAntStride +
                        polIdx * weightsPolStride +
                        frqIdx * weightsChnStride;

                    phasorAntennaCalibrations[outputIdx] = antennaCalibrationsCpp[inputIdx];
                }
            }
        }
    }

    State.inputBuffer = Duet<ArrayTensor<Device::CUDA, CI8>>(State.inputShape);

    State.pipelineRunner = std::make_shared<Runner>();
    BladePipeline::Config config = {
        .inputShape = State.inputShape,
        .outputShape = State.outputShape,

        .preBeamformerChannelizerRate = ata_a_config.channelizerRate,
        // .preBeamformerPolarizerConvertToCircular = BLADE_ATA_MODE_A_CIRCULAR_POLARIZATION,

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

        .beamformerIncoherentBeam = BLADE_ATA_MODE_A_OUTPUT_INCOHERENT_BEAM,

        .detectorEnable = true,
        .detectorIntegrationSize = ata_a_config.integrationSize,
        .detectorNumberOfOutputPolarizations = ata_a_config.numberOfOutputPolarizations,

        .castBlockSize = ata_a_config.castBlockSize,
        .channelizerBlockSize = ata_a_config.channelizerBlockSize,
        .beamformerBlockSize = ata_a_config.beamformerBlockSize,
        .detectorBlockSize = ata_a_config.detectorBlockSize
    };
    State.pipelineRunner->connect(
        State.pipeline,
        config,
        {
            .dut = blockDut1,
            .julianDate = blockJulianDate,
            .buffer = State.inputBuffer,
        }
    );

    State.InputPointerMap.reserve(State.pipelineRunner->numberOfStreams());
    State.OutputPointerMap.reserve(State.pipelineRunner->numberOfStreams());

    return true;
}

void blade_ata_a_terminate() {}

size_t blade_ata_a_get_input_size() {
    assert(State.pipelineRunner);
    return State.pipeline->getInputBuffer().shape().size();
}

size_t blade_ata_a_get_output_size() {
    assert(State.pipelineRunner);
    return State.pipeline->getOutputBuffer().shape().size();
}

void blade_ata_a_set_block_time_mjd(double mjd) {
    blockJulianDate.data()[0] = mjd;
}

void blade_ata_a_set_block_dut1(double dut1) {
    blockDut1.data()[0] = dut1;
}

void blade_ata_a_register_user_data(void* user_data) {
    State.UserData = user_data;
}

void blade_ata_a_register_input_buffer_prefetch_cb(blade_stateful_cb* f) {
    State.Callbacks.InputBufferPrefetch = f;
}

void blade_ata_a_register_input_buffer_fetch_cb(blade_input_buffer_fetch_cb* f) {
    State.Callbacks.InputBufferFetch = f;
}

void blade_ata_a_register_input_buffer_enqueued_cb(blade_input_buffer_enqueued_cb* f) {
    State.Callbacks.InputBufferEnqueued = f;
}

void blade_ata_a_register_input_buffer_ready_cb(blade_input_buffer_ready_cb* f) {
    State.Callbacks.InputBufferReady = f;
}

void blade_ata_a_register_output_buffer_fetch_cb(blade_output_buffer_fetch_cb* f) {
    State.Callbacks.OutputBufferFetch = f;
}

void blade_ata_a_register_output_buffer_ready_cb(blade_output_buffer_ready_cb* f) {
    State.Callbacks.OutputBufferReady = f;
}

void blade_ata_a_register_blade_queued_input_clear_cb(blade_clear_queued_cb* f) {
    State.Callbacks.InputClear = f;
}

void blade_ata_a_register_blade_queued_output_clear_cb(blade_clear_queued_cb* f) {
    State.Callbacks.OutputClear = f;
}

bool blade_ata_a_compute_step() {
    bool prefetch = State.Callbacks.InputBufferPrefetch(State.UserData);
    
    if(!State.pipelineRunner) {
        return false;
    }

    if (prefetch) {
        size_t bufferId_input;
        size_t bufferId_output;
        void* externalBuffer_input = nullptr;
        void* externalBuffer_output = nullptr;
        // Calls client callback to request empty input buffer.
        if (!State.Callbacks.InputBufferFetch(State.UserData, &externalBuffer_input, &bufferId_input)) {
            return false;
        }

        if (!State.Callbacks.OutputBufferFetch(State.UserData, &externalBuffer_output, &bufferId_output)) {
            BL_WARN("No output buffer available. Skipping input buffer {}.", bufferId_input);
            State.Callbacks.InputClear(State.UserData, bufferId_input);
            return false;
        }

        // Keeps track of pointer for "ready" callback.
        State.InputPointerMap.insert({bufferId_input, externalBuffer_input});
        State.OutputPointerMap.insert({bufferId_output, externalBuffer_output});

        // Create Memory::ArrayTensor from RAW pointer.
        auto input = ArrayTensor<Device::CPU, CI8>(externalBuffer_input, State.inputShape);
        auto output = ArrayTensor<Device::CPU, BLADE_ATA_MODE_A_OUTPUT_ELEMENT_T>(externalBuffer_output, State.pipeline->getOutputBuffer().shape());

        // Transfer input memory to the pipeline.
        auto inputCallback = [&](){
            BL_CHECK(State.pipelineRunner->copy(State.inputBuffer, input));
            return Result::SUCCESS;
        };
        auto outputCallback = [&](){
            BL_CHECK(State.pipelineRunner->copy(output, State.pipeline->getOutputBuffer()));
            return Result::SUCCESS;
        };
        State.pipelineRunner->enqueue(inputCallback, outputCallback, bufferId_input, bufferId_output);

        // Asynchronous CPU work
        State.Callbacks.InputBufferEnqueued(State.UserData, bufferId_input, bufferId_output);

        // Increment counter.
        State.StepCount++; 
    }

    // Dequeue last runner job and recycle output buffer.
    State.pipelineRunner->dequeue(
        [&](
            const U64& inputId, 
            const U64& outputId,
            const bool& didOutput
        ){
            if (didOutput) {
                const auto& recycleBuffer_input = State.OutputPointerMap[inputId];
                const auto& recycleBuffer_output = State.OutputPointerMap[outputId];
                
                State.Callbacks.InputBufferReady(State.UserData, recycleBuffer_input, inputId);
                State.Callbacks.OutputBufferReady(State.UserData, recycleBuffer_output, outputId);
                State.InputPointerMap.erase(inputId);
                State.OutputPointerMap.erase(outputId);
            }
            return Result::SUCCESS;
        }
    );

    // Return buffer was queued.
    return true;
}
