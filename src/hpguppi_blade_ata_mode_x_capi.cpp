#include <cassert>
#include <memory>

#include "blade/base.hh"
#include "blade/runner.hh"
#include "blade/bundles/generic/mode_x.hh"

extern "C" {
#include "hpguppi_blade_ata_mode_x_capi.h"
}

using namespace Blade;

using ModeX = Bundles::Generic::ModeX<CI8, BLADE_ATA_MODE_X_OUTPUT_ELEMENT_T>;

template<typename IT, typename OT>
class ModeXRunner : public Runner {
 public:
    struct Config {
        ArrayShape inputShape;
        ArrayShape outputShape;
    };

    explicit ModeXRunner(const Config& config, U64 channelizationRate, U64 integrationRate, U64 frequencyIntegrationRate)
        : inputBuffer(config.inputShape),
          outputBuffer(config.outputShape)
    {
        if (channelizationRate > 1 ) {
            // Spectral line mode (16 MHz firmware)
            if (channelizationRate < config.inputShape.numberOfTimeSamples()) {
                BL_FATAL("Channelizer rate must be a multiple of NTIME: {} < {}.", channelizationRate, config.inputShape.numberOfTimeSamples());
                throw Result::ASSERTION_ERROR;
            }
            if (channelizationRate%config.inputShape.numberOfTimeSamples() != 0) {
                BL_FATAL("Channelizer rate must be a multiple of NTIME: {}%{} != 0.", channelizationRate, config.inputShape.numberOfTimeSamples());
                throw Result::ASSERTION_ERROR;
            }
        }

        bool contiuumNotSpectral = channelizationRate == 1;
        // modeX will always produce T=1 output shapes:
        //   if channelizer is not bypassed, modeX will upchannelize all time first
        //   if channelizer is bypassed, modeX will integrate all input time
        // further integration happens thereafter to reach integrationRate
        U64 preCorrelatorStackerRate = contiuumNotSpectral ? 1 : channelizationRate/config.inputShape.numberOfTimeSamples();
        U64 correlatorIntegrationRate = contiuumNotSpectral ? integrationRate/config.inputShape.numberOfTimeSamples() : integrationRate;
        if (correlatorIntegrationRate % BLADE_ATA_MODE_X_INTEGRATION_FACTOR != 0) {
            BL_FATAL("Correlator integration rate must be a multiple of INTEGRATION_FACTOR: {}%{} != 0.", correlatorIntegrationRate, BLADE_ATA_MODE_X_INTEGRATION_FACTOR);
            throw Result::ASSERTION_ERROR;
        }

        ModeX::Config cfg = {
            .inputShape = config.inputShape,
            .outputShape = config.outputShape,
            .preChannelizerStackerMultiplier = 1,
            .channelizerBypass = contiuumNotSpectral,
            
            .preCorrelatorStackerMultiplier = BLADE_ATA_MODE_X_INTEGRATION_FACTOR,
            .correlatorIntegrationRate = correlatorIntegrationRate/BLADE_ATA_MODE_X_INTEGRATION_FACTOR,
            .correlatorConjugateAntennaIndex = BLADE_ATA_MODE_X_CONJUGATION_INDEX,

            .correlatorUseSharedMemory = false, //contiuumNotSpectral,
            .correlatorCalculationMode = CALC_MODE::INTEGER, // contiuumNotSpectral ? CALC_MODE::INTEGER : CALC_MODE::DOUBLE_PRECISION_FP,
            
            .postCorrelatorFrequencyIntegrationRate = frequencyIntegrationRate, // TODO support repeated integrations...

            .correlatorBlockSize = contiuumNotSpectral ? (U64) 64 : (U64) 32
        };
        this->connect(
            pipeline,
            cfg,
            {
                .buffer = inputBuffer
            }
        );
        this->compile();
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

 private:
    std::shared_ptr<ModeX> pipeline;

    Duet<ArrayTensor<Device::CUDA, IT>> inputBuffer;
    Duet<ArrayTensor<Device::CUDA, OT>> outputBuffer;
};

static struct {
    U64 enqueueCount = 0;
    void* UserData = nullptr;
    U64 InputProcessingRunningAverage;
    std::unordered_map<U64, struct timespec*> InputTimestampMap;
    std::unordered_map<U64, void*> InputPointerMap;
    std::unordered_map<U64, void*> OutputPointerMap;

    std::shared_ptr<ModeXRunner<CI8, BLADE_ATA_MODE_X_OUTPUT_ELEMENT_T>> pipelineRunner;

    size_t bufferId_output = 0;
    
    ArrayTensor<Device::CPU, CI8> debugInput;

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

bool blade_ata_x_initialize(
    struct blade_ata_mode_x_config ata_x_config,
    size_t numberOfWorkers,
    struct blade_ata_observation_meta* observationMeta,
    struct LonLatAlt* arrayReferencePosition,
    double* obs_phase_center_radecrad,
    double* beamCoordinates_radecrad,
    double* antennaPositions_xyz,
    double _Complex* antennaCalibrations
) {
    using namespace std::complex_literals;
    using ModeXRunner = ModeXRunner<CI8, BLADE_ATA_MODE_X_OUTPUT_ELEMENT_T>;

    BL_INFO("Initializing...");
    if (State.pipelineRunner) {
        BL_FATAL("Can't initialize because Blade Runner is already initialized.");
        throw Result::ASSERTION_ERROR;
    }

    State.inputShape = ArrayShape({
        ata_x_config.inputDims.NANTS,
        ata_x_config.inputDims.NCHANS,
        ata_x_config.inputDims.NTIME,
        ata_x_config.inputDims.NPOLS,
    });

    State.outputShape = ArrayShape({
        ata_x_config.inputDims.NANTS*(ata_x_config.inputDims.NANTS+1)/2,
        ata_x_config.inputDims.NCHANS*ata_x_config.channelizerRate/ata_x_config.frequencyIntegrationSize,
        1,
        ata_x_config.inputDims.NPOLS * ata_x_config.inputDims.NPOLS,
    });

    ModeXRunner::Config config = {
        .inputShape = State.inputShape,
        .outputShape = State.outputShape,
    };
    State.pipelineRunner = std::make_shared<ModeXRunner>(
        config,
        ata_x_config.channelizerRate,
        ata_x_config.integrationSize,
        ata_x_config.frequencyIntegrationSize
    );

    State.enqueueCount = 0;
    State.InputProcessingRunningAverage = 0;
    State.InputTimestampMap.reserve(State.pipelineRunner->numberOfStreams());
    State.InputPointerMap.reserve(State.pipelineRunner->numberOfStreams());
    State.OutputPointerMap.reserve(State.pipelineRunner->numberOfStreams());

    // State.debugInput = ArrayTensor<Device::CPU, CI8>(State.inputShape);
    // size_t index = 0;
    // for (int a = 0; a < ata_x_config.inputDims.NANTS; a++) {
    //     for (int c = 0; c < ata_x_config.inputDims.NCHANS; c++) {
    //         for (int t = 0; t < ata_x_config.inputDims.NTIME; t++) {
    //             for (int p = 0; p < ata_x_config.inputDims.NPOLS; p++) {
    //                 if (a == 2) {
    //                     if (p==0)
    //                         State.debugInput[index++] = std::complex<int8_t>(1, 0);
    //                     else
    //                         State.debugInput[index++] = std::complex<int8_t>(0, 0);
    //                 } else if (a == 3) {
    //                     if (p==0)
    //                         State.debugInput[index++] = std::complex<int8_t>(0, 0);
    //                     else
    //                         State.debugInput[index++] = std::complex<int8_t>(1, 0);
    //                 } else if (a == 1) {
    //                     State.debugInput[index++] = std::complex<int8_t>(1, 0);
    //                 } else {
    //                     State.debugInput[index++] = std::complex<int8_t>(a*3, (t%250-125)+p+1);
    //                 }
    //             }
    //         }
    //     }
    // }
    return true;
}

void blade_ata_x_terminate() {
    if (!State.pipelineRunner) {
        BL_WARN("No pipeline to terminate.")
        return;
    }
    State.pipelineRunner.reset();
    if (State.pipelineRunner) {
        BL_WARN("Reset didn't invalidate pipelineRunner.")
    }
    // State.pipeline.reset();
    State.InputTimestampMap.clear();
    for (const auto& [inputId, recycleBuffer_input] : State.InputPointerMap) {
        State.Callbacks.InputBufferReady(State.UserData, recycleBuffer_input, inputId);
    }
    State.InputPointerMap.clear();
    for (const auto& [outputId, recycleBuffer_output] : State.OutputPointerMap) {
        State.Callbacks.OutputBufferReady(State.UserData, recycleBuffer_output, outputId);
    }
    State.OutputPointerMap.clear();
}

size_t blade_ata_x_get_input_size() {
    assert(State.pipelineRunner);
    return State.inputShape.size();
}

size_t blade_ata_x_get_output_size() {
    assert(State.pipelineRunner);
    return State.outputShape.size();
}

void blade_ata_x_set_block_time_mjd(double mjd) {
    // blockJulianDate.data()[0] = mjd;
}

void blade_ata_x_set_block_dut1(double dut1) {
    // blockDut1.data()[0] = dut1;
}

void blade_ata_x_register_user_data(void* user_data) {
    State.UserData = user_data;
}

void blade_ata_x_register_input_buffer_prefetch_cb(blade_stateful_cb* f) {
    State.Callbacks.InputBufferPrefetch = f;
}

void blade_ata_x_register_input_buffer_fetch_cb(blade_input_buffer_fetch_cb* f) {
    State.Callbacks.InputBufferFetch = f;
}

void blade_ata_x_register_input_buffer_enqueued_cb(blade_input_buffer_enqueued_cb* f) {
    State.Callbacks.InputBufferEnqueued = f;
}

void blade_ata_x_register_input_buffer_ready_cb(blade_input_buffer_ready_cb* f) {
    State.Callbacks.InputBufferReady = f;
}

void blade_ata_x_register_output_buffer_fetch_cb(blade_output_buffer_fetch_cb* f) {
    State.Callbacks.OutputBufferFetch = f;
}

void blade_ata_x_register_output_buffer_ready_cb(blade_output_buffer_ready_cb* f) {
    State.Callbacks.OutputBufferReady = f;
}

void blade_ata_x_register_blade_queued_input_clear_cb(blade_clear_queued_cb* f) {
    State.Callbacks.InputClear = f;
}

void blade_ata_x_register_blade_queued_output_clear_cb(blade_clear_queued_cb* f) {
    State.Callbacks.OutputClear = f;
}

bool blade_ata_x_compute_step() {
    bool prefetch = State.Callbacks.InputBufferPrefetch(State.UserData);
    
    if(!State.pipelineRunner) {
        if (prefetch) {
            BL_WARN("Prefetch is true but pipelineRunner is null...");
        }
        return false;
    }

    if (prefetch) {
        size_t bufferId_input;
        void* externalBuffer_input = nullptr;
        // Calls client callback to request empty input buffer.
        if (!State.Callbacks.InputBufferFetch(State.UserData, &externalBuffer_input, &bufferId_input)) {
            BL_WARN("No input buffer retrieved.");
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
        auto output = ArrayTensor<Device::CPU, BLADE_ATA_MODE_X_OUTPUT_ELEMENT_T>(State.OutputPointerMap[State.bufferId_output], State.outputShape);

        // Transfer input memory to the pipeline.
        auto inputCallback = [&](){
            // return State.pipelineRunner->transferIn(State.debugInput);
            return State.pipelineRunner->transferIn(input);
        };
        auto transferCallback = [&](){
            void* recycleBuffer_input = State.InputPointerMap[bufferId_input];
            State.Callbacks.InputBufferReady(State.UserData, recycleBuffer_input, bufferId_input);
            State.InputPointerMap.erase(bufferId_input);
            return Result::SUCCESS;
        };
        auto resultCallback = [&](){
            return State.pipelineRunner->transferResult();
        };
        auto outputCallback = [&](){
            return State.pipelineRunner->transferOut(output);
        };

        // Dequeue last runner job and recycle output buffer.
        // only blocks if queue is full which would cause enqueue failure anyway...
        State.pipelineRunner->dequeue(
            [&](
                const U64& inputId, 
                const U64& outputId,
                const bool& didOutput
            ){
                // void* recycleBuffer_input = State.InputPointerMap[inputId];
                // State.Callbacks.InputBufferReady(State.UserData, recycleBuffer_input, inputId);
                // State.InputPointerMap.erase(inputId);

                if (didOutput) { // should assert this really
                    void* recycleBuffer_output = State.OutputPointerMap[outputId];
                    State.Callbacks.OutputBufferReady(State.UserData, recycleBuffer_output, outputId);
                    State.OutputPointerMap.erase(outputId);
                }
                return Result::SUCCESS;
            }
        );

        if ( Result::SUCCESS !=
            State.pipelineRunner->enqueue(inputCallback, transferCallback, resultCallback, outputCallback, bufferId_input, State.bufferId_output)
        ) {
            BL_FATAL("Could not enqueue block #{}!", bufferId_input);
        }
        // Asynchronous CPU work
        State.Callbacks.InputBufferEnqueued(State.UserData, bufferId_input, State.bufferId_output);        
    }


    // Return buffer was queued.
    return true;
}
