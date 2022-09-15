#ifndef HPGUPPI_BLADE
#define HPGUPPI_BLADE

typedef bool (blade_stateful_cb)(void*);
typedef bool (blade_input_buffer_fetch_cb)(void*, void**, size_t*);
typedef void (blade_input_buffer_enqueued_cb)(void*, size_t);
typedef void (blade_input_buffer_ready_cb)(void*, const void*, size_t);
typedef bool (blade_output_buffer_fetch_cb)(void*, void**, size_t*);
typedef void (blade_output_buffer_ready_cb)(void*, const void*, size_t);
typedef void (blade_clear_queued_cb)(void*, size_t);

bool blade_use_device(int device_id);
bool blade_pin_memory(void* buffer, size_t size);

#endif