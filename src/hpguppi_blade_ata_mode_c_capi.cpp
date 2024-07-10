#include <cassert>
#include <memory>

#include "blade/base.hh"
#include "blade/runner.hh"
#include "blade/modules/correlator.hh"
#include "blade/modules/cast.hh"

extern "C" {
#include "hpguppi_blade_ata_mode_c_capi.h"
}

using namespace Blade;

using Cast = Modules::Cast<CI8, CF32>;
using Correlator = Modules::Correlator<CF32, CF32>;

static struct {
    U64 StepCount = 0;
    void* UserData = nullptr;
    std::unordered_map<U64, void*> InputPointerMap;
    std::unordered_map<U64, void*> OutputPointerMap;

    std::shared_ptr<Runner> pipelineRunner;
    std::shared_ptr<Cast> inputCast;
    std::shared_ptr<Correlator> correlator;
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

bool blade_ata_c_initialize(
    struct blade_ata_mode_c_config ata_c_config,
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
    BL_INFO("Initializing...");

    State.inputShape = ArrayShape({
        ata_c_config.inputDims.NANTS,
        ata_c_config.inputDims.NCHANS,
        ata_c_config.inputDims.NTIME,
        ata_c_config.inputDims.NPOLS,
    });

    auto outputShape = ArrayShape({
        ata_c_config.inputDims.NANTS*(ata_c_config.inputDims.NANTS+1)/2,
        ata_c_config.inputDims.NCHANS,
        ata_c_config.inputDims.NTIME / ata_c_config.integrationSize,
        ata_c_config.inputDims.NPOLS * ata_c_config.inputDims.NPOLS,
    });

    State.inputBuffer = Duet<ArrayTensor<Device::CUDA, CI8>>(State.inputShape);

    State.pipelineRunner = std::make_shared<Runner>();
    Cast::Config config = {
        .blockSize = ata_c_config.castBlockSize,
    };
    State.pipelineRunner->connect(
        State.inputCast,
        config,
        {
            .buf = State.inputBuffer,
        }
    );

    Correlator::Config configCorrelator = {
        .integrationSize = ata_c_config.integrationSize,
        .blockSize = ata_c_config.correlatorBlockSize,
    };
    State.pipelineRunner->connect(
        State.correlator,
        configCorrelator,
        {
            .buf = State.inputCast->getOutputBuffer(),
        }
    );

    State.InputPointerMap.reserve(State.pipelineRunner->numberOfStreams());
    State.OutputPointerMap.reserve(State.pipelineRunner->numberOfStreams());

    return true;
}

void blade_ata_c_terminate() {}

size_t blade_ata_c_get_input_size() {
    assert(State.inputCast);
    return State.inputCast->getInputBuffer().shape().size();
}

size_t blade_ata_c_get_output_size() {
    assert(State.correlator);
    return State.correlator->getOutputBuffer().shape().size();
}

void blade_ata_c_set_block_time_mjd(double mjd) {
    // blockJulianDate.data()[0] = mjd;
}

void blade_ata_c_set_block_dut1(double dut1) {
    // blockDut1.data()[0] = dut1;
}

void blade_ata_c_register_user_data(void* user_data) {
    State.UserData = user_data;
}

void blade_ata_c_register_input_buffer_prefetch_cb(blade_stateful_cb* f) {
    State.Callbacks.InputBufferPrefetch = f;
}

void blade_ata_c_register_input_buffer_fetch_cb(blade_input_buffer_fetch_cb* f) {
    State.Callbacks.InputBufferFetch = f;
}

void blade_ata_c_register_input_buffer_enqueued_cb(blade_input_buffer_enqueued_cb* f) {
    State.Callbacks.InputBufferEnqueued = f;
}

void blade_ata_c_register_input_buffer_ready_cb(blade_input_buffer_ready_cb* f) {
    State.Callbacks.InputBufferReady = f;
}

void blade_ata_c_register_output_buffer_fetch_cb(blade_output_buffer_fetch_cb* f) {
    State.Callbacks.OutputBufferFetch = f;
}

void blade_ata_c_register_output_buffer_ready_cb(blade_output_buffer_ready_cb* f) {
    State.Callbacks.OutputBufferReady = f;
}

void blade_ata_c_register_blade_queued_input_clear_cb(blade_clear_queued_cb* f) {
    State.Callbacks.InputClear = f;
}

void blade_ata_c_register_blade_queued_output_clear_cb(blade_clear_queued_cb* f) {
    State.Callbacks.OutputClear = f;
}

bool blade_ata_c_compute_step() {
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
        auto output = ArrayTensor<Device::CPU, CF32>(externalBuffer_output, State.correlator->getOutputBuffer().shape());

        // Transfer input memory to the pipeline.
        auto inputCallback = [&](){
            BL_CHECK(State.pipelineRunner->copy(State.inputBuffer, input));
            return Result::SUCCESS;
        };
        auto outputCallback = [&](){
            BL_CHECK(State.pipelineRunner->copy(output, State.correlator->getOutputBuffer()));
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
