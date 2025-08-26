#ifndef BLADE_ATA_MODE_K_H
#define BLADE_ATA_MODE_K_H

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <complex.h>

#include "hpguppi_blade_capi.h"
#include "hpguppi_blade_ata_mode_k_config.h"

struct blade_ata_mode_k_config {
    uint32_t channelizerRate;
    uint32_t beamformerBeams;
    uint32_t integrationSize;

    bool kurtosisEnabled;
    int32_t kurtosisSigma;
    uint32_t kurtosisChannelLength;
    uint32_t kurtosisNumberOfMaskRuns;
    char kurtosisMaskOutputFilepath[256];

    uint32_t castBlockSize;
};

static const struct blade_ata_mode_k_config BLADE_ATA_MODE_K_CONFIG = {
    1, // .channelizerRate
    1, // .beamformerBeams
    1, // .integrationSize

    false, // .kurtosisEnabled/
    3, // .kurtosisSigma
    256, // .kurtosisChannelLength
    128, // .kurtosisNumberOfMaskRuns
    "/dev/null", // .kurtosisMaskOutputFilepath

    512, // .castBlockSize
};

bool blade_ata_k_initialize(
    struct blade_ata_mode_k_config ata_k_config
);
size_t blade_ata_k_get_output_size();

void blade_ata_k_set_block_time_mjd(double mjd);
void blade_ata_k_set_block_dut1(double dut1);

void blade_ata_k_register_user_data(void* user_data);
void blade_ata_k_register_input_buffer_prefetch_cb(blade_stateful_cb* f);
void blade_ata_k_register_input_buffer_fetch_cb(blade_input_buffer_fetch_cb* f);
void blade_ata_k_register_input_buffer_enqueued_cb(blade_input_buffer_enqueued_cb* f);
void blade_ata_k_register_input_buffer_ready_cb(blade_input_buffer_ready_cb* f);
void blade_ata_k_register_output_buffer_fetch_cb(blade_output_buffer_fetch_cb* f);
void blade_ata_k_register_output_buffer_ready_cb(blade_output_buffer_ready_cb* f);

void blade_ata_k_register_blade_queued_input_clear_cb(blade_clear_queued_cb* f);
void blade_ata_k_register_blade_queued_output_clear_cb(blade_clear_queued_cb* f);

bool blade_ata_k_compute_step();

void blade_ata_k_terminate();

#endif
