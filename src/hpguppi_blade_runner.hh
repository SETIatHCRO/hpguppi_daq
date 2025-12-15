#ifndef HPGUPPI_BLADE_RUNNER
#define HPGUPPI_BLADE_RUNNER

#include "hpguppi_blade_capi.h"

#include "blade/base.hh"
#include "blade/runner.hh"

using namespace Blade;

class ModeRunner : public Runner {
    public:
        virtual Result transferIn(const ArrayTensor<Device::CPU, CI8>& cpuInputBuffer) = 0;
        virtual Result transferResult() = 0;
        virtual Result transferOut(void* cpuOutputBuffer) = 0;
        virtual size_t outputByteSize() = 0;
};

typedef struct {
    U64 StepCount = 0;
    void* UserData = nullptr;
    U64 InputProcessingRunningAverage;
    std::unordered_map<U64, struct timespec*> InputTimestampMap;
    std::unordered_map<U64, void*> InputPointerMap;
    std::unordered_map<U64, void*> OutputPointerMap;
    
    std::shared_ptr<ModeRunner> pipelineRunner;
    
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
        blade_reset_output_index_cb* OutputBufferResetIndex;
    } Callbacks;
} StateStruct;

extern StateStruct State;

#endif