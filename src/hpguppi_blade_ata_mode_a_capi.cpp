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

template<typename IT, typename OT>
class ModeARunner : public Runner {
 public:
    struct Config {
        ArrayShape inputShape;
        ArrayShape outputShape;
    };

    explicit ModeARunner(
        const Config& config,
        struct blade_ata_mode_a_config* ata_a_config,
        struct blade_ata_observation_meta* observationMeta,
        struct LonLatAlt* arrayReferencePosition,
        double* obs_phase_center_radecrad,
        double* beamCoordinates_radecrad,
        double* antennaPositions_xyz,
        double _Complex* antennaCalibrations
    )
        : inputBuffer(config.inputShape),
          outputBuffer(config.outputShape),
          blockJulianDate({1}),
          blockDut1({1})
    {
        std::vector<XYZ> antennaPositions(ata_a_config->inputDims.NANTS);
        std::vector<RA_DEC> beamCoordinates(ata_a_config->beamformerBeams);
        std::vector<std::complex<double>> antennaCalibrationsCpp(
                ata_a_config->inputDims.NANTS*\
                ata_a_config->inputDims.NCHANS*\
                ata_a_config->inputDims.NPOLS);
        int i;
        for(i = 0; i < ata_a_config->inputDims.NANTS; i++){
            antennaPositions[i].X = antennaPositions_xyz[i*3 + 0];
            antennaPositions[i].Y = antennaPositions_xyz[i*3 + 1];
            antennaPositions[i].Z = antennaPositions_xyz[i*3 + 2];
        }
        for(i = 0; i < ata_a_config->beamformerBeams; i++){
            beamCoordinates[i].RA = beamCoordinates_radecrad[i*2 + 0];
            beamCoordinates[i].DEC = beamCoordinates_radecrad[i*2 + 1];
        }
        memcpy(antennaCalibrationsCpp.data(), antennaCalibrations,
                antennaCalibrationsCpp.size()*sizeof(antennaCalibrationsCpp[0]));

        auto phasorAntennaCalibrations = ArrayTensor<Device::CPU, CF64>({
            ata_a_config->inputDims.NANTS,
            ata_a_config->inputDims.NCHANS * ata_a_config->channelizerRate,
            1,
            ata_a_config->inputDims.NPOLS,
        });

        const size_t calAntStride = 1;
        const size_t calPolStride = config.inputShape.numberOfAspects() * calAntStride;
        const size_t calChnStride = config.inputShape.numberOfPolarizations() * calPolStride;

        const size_t weightsPolStride = 1;
        const size_t weightsChnStride = config.inputShape.numberOfPolarizations() * weightsPolStride;
        const size_t weightsAntStride = config.inputShape.numberOfFrequencyChannels() * ata_a_config->channelizerRate * weightsChnStride;
        BL_INFO("Expanding the {} coarse-channel coefficients by a factor of {}.", config.inputShape.numberOfFrequencyChannels(), ata_a_config->channelizerRate);

        U64 inputIdx, frqIdx, outputIdx, antIdx, chnIdx, polIdx, fchIdx;
        for (antIdx = 0; antIdx < config.inputShape.numberOfAspects(); antIdx++) {
            for (chnIdx = 0; chnIdx < config.inputShape.numberOfFrequencyChannels(); chnIdx++) {
                for (polIdx = 0; polIdx < config.inputShape.numberOfPolarizations(); polIdx++) {
                    inputIdx = chnIdx * calChnStride +
                        polIdx * calPolStride + 
                        antIdx * calAntStride;
                    for (fchIdx = 0; fchIdx < ata_a_config->channelizerRate; fchIdx++) {
                        frqIdx = chnIdx * ata_a_config->channelizerRate + fchIdx;
                        outputIdx = antIdx * weightsAntStride +
                            polIdx * weightsPolStride +
                            frqIdx * weightsChnStride;

                        phasorAntennaCalibrations[outputIdx] = antennaCalibrationsCpp[inputIdx];
                    }
                }
            }
        }

        this->inputBuffer = Duet<ArrayTensor<Device::CUDA, CI8>>(config.inputShape);
        BladePipeline::Config pipeline_config = {
            .inputShape = config.inputShape,
            .outputShape = config.outputShape,

            .preBeamformerChannelizerRate = ata_a_config->channelizerRate,
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
            .detectorIntegrationRate = ata_a_config->integrationSize,
            .detectorNumberOfOutputPolarizations = ata_a_config->numberOfOutputPolarizations,

            .casterBlockSize = ata_a_config->castBlockSize,
            .channelizerBlockSize = ata_a_config->channelizerBlockSize,
            .beamformerBlockSize = ata_a_config->beamformerBlockSize,
            .detectorBlockSize = ata_a_config->detectorBlockSize
        };
        this->connect(
            this->pipeline,
            pipeline_config,
            {
                .dut = blockDut1,
                .julianDate = blockJulianDate,
                .buffer = this->inputBuffer,
            }
        );
    }


    Result transferIn(const ArrayTensor<Device::CPU, IT>& cpuInputBuffer) {
        BL_CHECK(this->copy(inputBuffer, cpuInputBuffer));
        return Result::SUCCESS;
    }

    Result transferResult() {
        BL_CHECK(this->copy(outputBuffer, pipeline->getOutputBuffer()));
        return Result::SUCCESS;
    }
    Result transferOut(ArrayTensor<Device::CPU, OT>& cpuOutputBuffer) {
        BL_CHECK(this->copy(cpuOutputBuffer, outputBuffer));
        return Result::SUCCESS;
    }

    void setJulianDate(F64 value) {
        this->blockJulianDate.data()[0] = value;
    }

    void setDut1(F64 value) {
        this->blockDut1.data()[0] = value;
    }

 private:
    std::shared_ptr<Runner> pipelineRunner;
    std::shared_ptr<BladePipeline> pipeline;
    
    Duet<ArrayTensor<Device::CUDA, IT>> inputBuffer;
    Duet<ArrayTensor<Device::CUDA, OT>> outputBuffer;
    Tensor<Device::CPU, F64> blockJulianDate;
    Tensor<Device::CPU, F64> blockDut1;
};

static struct {
    U64 StepCount = 0;
    void* UserData = nullptr;
    std::unordered_map<U64, void*> InputPointerMap;
    std::unordered_map<U64, void*> OutputPointerMap;
    
    std::shared_ptr<ModeARunner<CI8, BLADE_ATA_MODE_A_OUTPUT_ELEMENT_T>> pipelineRunner;
    
    size_t bufferId_output = 0;

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
    using ModeARunner = ModeARunner<CI8, BLADE_ATA_MODE_A_OUTPUT_ELEMENT_T>;

    if (State.pipelineRunner) {
        BL_FATAL("Can't initialize because Blade Runner is already initialized.");
        throw Result::ASSERTION_ERROR;
    }
    BL_INFO("Initializing...");

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

    ModeARunner::Config config = {
        .inputShape = State.inputShape,
        .outputShape = State.outputShape
    };
    State.pipelineRunner = std::make_shared<ModeARunner>(
        config,
        &ata_a_config,
        observationMeta,
        arrayReferencePosition,
        obs_phase_center_radecrad,
        beamCoordinates_radecrad,
        antennaPositions_xyz,
        antennaCalibrations
    );
    State.InputPointerMap.reserve(State.pipelineRunner->numberOfStreams());
    State.OutputPointerMap.reserve(State.pipelineRunner->numberOfStreams());

    return true;
}

void blade_ata_a_terminate() {
    if (!State.pipelineRunner) {
        BL_WARN("No pipeline to terminate.")
        return;
    }
    State.pipelineRunner.reset();

    for (const auto& [inputId, recycleBuffer_input] : State.InputPointerMap) {
        State.Callbacks.InputClear(State.UserData, inputId);
    }
    State.InputPointerMap.clear();
    for (const auto& [outputId, recycleBuffer_output] : State.OutputPointerMap) {
        State.Callbacks.OutputClear(State.UserData, outputId);
    }
    State.OutputPointerMap.clear();
}

size_t blade_ata_a_get_input_size() {
    assert(State.pipelineRunner);
    return State.inputShape.size();
}

size_t blade_ata_a_get_output_size() {
    assert(State.pipelineRunner);
    return State.outputShape.size();
}

void blade_ata_a_set_block_time_mjd(double mjd) {
    State.pipelineRunner->setJulianDate(mjd);
}

void blade_ata_a_set_block_dut1(double dut1) {
    State.pipelineRunner->setDut1(dut1);
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
        void* externalBuffer_input = nullptr;
        // Calls client callback to request empty input buffer.
        if (!State.Callbacks.InputBufferFetch(State.UserData, &externalBuffer_input, &bufferId_input)) {
            return false;
        }

        if (State.pipelineRunner->computeCurrentStepCount() == 0) {
            void* externalBuffer_output = nullptr;
            if (!State.Callbacks.OutputBufferFetch(State.UserData, &externalBuffer_output, &State.bufferId_output)) {
                BL_WARN("No output buffer available. Skipping input buffer {}.", bufferId_input);
                State.Callbacks.InputClear(State.UserData, bufferId_input);
                return false;
            }
            State.OutputPointerMap.insert({State.bufferId_output, externalBuffer_output});
        }

        // Create Memory::ArrayTensor from RAW pointer.
        auto input = ArrayTensor<Device::CPU, CI8>(externalBuffer_input, State.inputShape);
        State.InputPointerMap.insert({bufferId_input, externalBuffer_input});
        auto output = ArrayTensor<Device::CPU, BLADE_ATA_MODE_A_OUTPUT_ELEMENT_T>(State.OutputPointerMap[State.bufferId_output], State.outputShape);

        // Transfer input memory to the pipeline.
        auto inputCallback = [&](){
            return State.pipelineRunner->transferIn(input);
        };
        auto resultCallback = [&](){
            return State.pipelineRunner->transferResult();
        };
        auto outputCallback = [&](){
            return State.pipelineRunner->transferOut(output);
        };

        if ( Result::SUCCESS !=
            State.pipelineRunner->enqueue(inputCallback, resultCallback, outputCallback, bufferId_input, State.bufferId_output)
        ) {
            // Dequeue last runner job and recycle output buffer.
            State.pipelineRunner->dequeue(
                [&](
                    const U64& inputId, 
                    const U64& outputId,
                    const bool& didOutput
                ){
                    void* recycleBuffer_input = State.InputPointerMap[inputId];
                    State.Callbacks.InputBufferReady(State.UserData, recycleBuffer_input, inputId);
                    State.InputPointerMap.erase(inputId);

                    if (didOutput) { // should assert this really
                    
                        void* recycleBuffer_output = State.OutputPointerMap[outputId];
                        State.Callbacks.OutputBufferReady(State.UserData, recycleBuffer_output, outputId);
                        State.OutputPointerMap.erase(outputId);
                    }
                    return Result::SUCCESS;
                }
            );
            State.pipelineRunner->enqueue(inputCallback, resultCallback, outputCallback, bufferId_input, State.bufferId_output);// should assert or something
        }
        // Asynchronous CPU work
        State.Callbacks.InputBufferEnqueued(State.UserData, bufferId_input, State.bufferId_output);
    }
    // else {
    //     // Dequeue last runner job and recycle output buffer.
    //     State.pipelineRunner->dequeue(
    //         [&](
    //             const U64& inputId, 
    //             const U64& outputId,
    //             const bool& didOutput
    //         ){
    //             if (didOutput) {
    //                 void* recycleBuffer_output = State.OutputPointerMap[outputId];
    //                 State.Callbacks.OutputBufferReady(State.UserData, recycleBuffer_output, outputId);
    //                 State.OutputPointerMap.erase(outputId);
    //                 return Result::SUCCESS;
    //             }
    //             return Result::ERROR;
    //         }
    //     );
    // }

    // Return buffer was queued.
    return true;
}
