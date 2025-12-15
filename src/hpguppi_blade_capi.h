#ifndef HPGUPPI_BLADE
#define HPGUPPI_BLADE

#include <stdbool.h>
#include <stdint.h>
#include "config_input_shape.h"
#include "hpguppi_blade_ata_structs.h"

typedef bool (blade_stateful_cb)(void*);
typedef bool (blade_input_buffer_fetch_cb)(void*, void**, size_t*);
typedef void (blade_input_buffer_enqueued_cb)(void*, size_t, size_t);
typedef void (blade_input_buffer_ready_cb)(void*, const void*, size_t);
typedef bool (blade_output_buffer_fetch_cb)(void*, void**, size_t*);
typedef void (blade_output_buffer_ready_cb)(void*, const void*, size_t);
typedef void (blade_clear_queued_cb)(void*, size_t);
typedef void (blade_reset_output_index_cb)(void*);

bool blade_use_device(int device_id);
bool blade_pin_memory(void* buffer, size_t size);

void blade_terminate();
void blade_set_input_dimensions(
    struct blade_ata_input_dims* inputDims
);
size_t blade_get_input_size();
size_t blade_get_output_size();
size_t blade_get_output_byte_size();
bool blade_compute_step();

void blade_ata_register_user_data(void* user_data);
void blade_ata_register_input_buffer_prefetch_cb(blade_stateful_cb* f);
void blade_ata_register_input_buffer_fetch_cb(blade_input_buffer_fetch_cb* f);
void blade_ata_register_input_buffer_enqueued_cb(blade_input_buffer_enqueued_cb* f);
void blade_ata_register_input_buffer_ready_cb(blade_input_buffer_ready_cb* f);
void blade_ata_register_output_buffer_fetch_cb(blade_output_buffer_fetch_cb* f);
void blade_ata_register_output_buffer_ready_cb(blade_output_buffer_ready_cb* f);
void blade_ata_register_blade_queued_input_clear_cb(blade_clear_queued_cb* f);
void blade_ata_register_blade_reset_output_index_cb(blade_reset_output_index_cb* f);

#endif