#include <cassert>
#include <memory>

#include "blade/base.hh"

extern "C" {
    #include "hpguppi_blade_capi.h"
}

#include "hpguppi_blade_runner.hh"

StateStruct State;

using namespace Blade;

bool blade_use_device(int device_id) {
    return SetCudaDevice(device_id) == Result::SUCCESS;
}

bool blade_pin_memory(void* buffer, size_t size) {
    // return PageLock(ArrayTensor<Device::CPU, I8>(buffer, {size})) == Result::SUCCESS;
    cudaError_t val = cudaHostRegister(buffer, size, cudaHostRegisterDefault);
    if (val != cudaSuccess) {
        const char* err = cudaGetErrorString(val);
        BL_WARN("Failed to register CPU memory: {}", err);
        return false;

    }
    return true;
}

void blade_terminate() {
    if (!State.pipelineRunner) {
        BL_WARN("No pipeline to terminate.")
        return;
    }
    State.pipelineRunner.reset();
    if (State.pipelineRunner) {
        BL_WARN("Reset didn't invalidate pipelineRunner.")
    }

    for (const auto& [inputId, recycleBuffer_input] : State.InputPointerMap) {
        State.Callbacks.InputClear(State.UserData, inputId);
    }
    State.InputPointerMap.clear();

    State.OutputPointerMap.clear();
    State.Callbacks.OutputBufferResetIndex(State.UserData);
}

void blade_set_input_dimensions(
    struct blade_ata_input_dims* inputDims
) {
    assert(!State.pipelineRunner);
    State.inputShape = ArrayShape({
        inputDims->NANTS,
        inputDims->NCHANS,
        inputDims->NTIME,
        inputDims->NPOLS
    });
}

size_t blade_get_input_size() {
    assert(State.pipelineRunner);
    return State.inputShape.size();
}

size_t blade_get_output_size() {
    assert(State.pipelineRunner);
    return State.outputShape.size();
}

size_t blade_get_output_byte_size() {
    assert(State.pipelineRunner);
    return State.pipelineRunner->outputByteSize();
}

bool blade_compute_step() {
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
        // auto output = ArrayTensor<Device::CPU, BLADE_ATA_MODE_X_OUTPUT_ELEMENT_T>(State.OutputPointerMap[State.bufferId_output], State.outputShape);

        // Transfer input memory to the pipeline.
        auto inputCallback = [&](){
            return State.pipelineRunner->transferIn(input);
        };
        auto transferCallback = [&](){
            struct timespec ts_input = {0};
            clock_gettime(CLOCK_MONOTONIC, &ts_input);
            State.InputTimestampMap[bufferId_input] = &ts_input;
    
            void* recycleBuffer_input = State.InputPointerMap[bufferId_input];
            State.Callbacks.InputBufferReady(State.UserData, recycleBuffer_input, bufferId_input);
            State.InputPointerMap.erase(bufferId_input);
            return Result::SUCCESS;
        };
        auto resultCallback = [&](){
            return State.pipelineRunner->transferResult();
        };
        auto outputCallback = [&](){
            return State.pipelineRunner->transferOut(State.OutputPointerMap[State.bufferId_output]);
        };

        // Dequeue last runner job and recycle output buffer.
        // only blocks if queue is full which would cause enqueue failure anyway...
        State.pipelineRunner->dequeue(
            [&](
                const U64& inputId, 
                const U64& outputId,
                const bool& didOutput
            ){
                struct timespec ts_output = {0};
                clock_gettime(CLOCK_MONOTONIC, &ts_output);
                struct timespec* ts_input = State.InputTimestampMap[inputId];
                
                // BL_INFO("Block turnaround time: {} ns", ((int64_t)(ts_output.tv_sec-ts_input->tv_sec))*1000*1000*1000+(ts_output.tv_nsec-ts_input->tv_nsec));
                State.InputTimestampMap.erase(inputId);
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

void blade_ata_register_user_data(void* user_data) {
    State.UserData = user_data;
}

void blade_ata_register_input_buffer_prefetch_cb(blade_stateful_cb* f) {
    State.Callbacks.InputBufferPrefetch = f;
}

void blade_ata_register_input_buffer_fetch_cb(blade_input_buffer_fetch_cb* f) {
    State.Callbacks.InputBufferFetch = f;
}

void blade_ata_register_input_buffer_enqueued_cb(blade_input_buffer_enqueued_cb* f) {
    State.Callbacks.InputBufferEnqueued = f;
}

void blade_ata_register_input_buffer_ready_cb(blade_input_buffer_ready_cb* f) {
    State.Callbacks.InputBufferReady = f;
}

void blade_ata_register_output_buffer_fetch_cb(blade_output_buffer_fetch_cb* f) {
    State.Callbacks.OutputBufferFetch = f;
}

void blade_ata_register_output_buffer_ready_cb(blade_output_buffer_ready_cb* f) {
    State.Callbacks.OutputBufferReady = f;
}

void blade_ata_register_blade_queued_input_clear_cb(blade_clear_queued_cb* f) {
    State.Callbacks.InputClear = f;
}

void blade_ata_register_blade_reset_output_index_cb(blade_reset_output_index_cb* f) {
    State.Callbacks.OutputBufferResetIndex = f;
}