#ifndef BLADE_ATA_MODE_C_H
#define BLADE_ATA_MODE_C_H

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <complex.h>

#include "hpguppi_blade_capi.h"
#include "hpguppi_blade_ata_structs.h"
#include "hpguppi_blade_ata_mode_c_config.h"

struct blade_ata_mode_c_config {
    struct blade_ata_input_dims inputDims;
    uint32_t channelizerRate;
    uint32_t beamformerBeams;
    uint32_t integrationSize;

    uint32_t outputMemWidth;
    uint32_t outputMemPad;

    uint32_t castBlockSize;
    uint32_t correlatorBlockSize;
};

static const struct blade_ata_mode_c_config BLADE_ATA_MODE_C_CONFIG = {
    { // .inputDims
        1, // .NBEAMS
        BLADE_ATA_MODE_C_INPUT_NANT, // .NANTS 
        BLADE_ATA_MODE_C_ANT_NCHAN, // .NCHANS
        BLADE_ATA_MODE_C_NTIME, // .NTIME 
        BLADE_ATA_MODE_C_NPOL, // .NPOLS 
    }, // .inputDims
    1, // .channelizerRate
    1, // .beamformerBeams
    BLADE_ATA_MODE_C_INTEGRATION_SIZE, // .integrationSize

    BLADE_ATA_MODE_C_OUTPUT_MEMCPY2D_WIDTH, // .outputMemWidth
    BLADE_ATA_MODE_C_OUTPUT_MEMCPY2D_PAD, // .outputMemPad

    512, // .castBlockSize
    512, // .correlatorBlockSize
};

bool blade_ata_c_initialize(
    struct blade_ata_mode_c_config ata_c_config,
    size_t numberOfWorkers,
    struct blade_ata_observation_meta* observationMeta,
    struct LonLatAlt* arrayReferencePosition,
    double* obs_phase_center_radecrad,
    double* beamCoordinates_radecrad,
    double* antennaPositions_xyz,
    double _Complex* antennaCalibrations
);
size_t blade_ata_c_get_output_size();

void blade_ata_c_set_block_time_mjd(double mjd);
void blade_ata_c_set_block_dut1(double dut1);

void blade_ata_c_register_user_data(void* user_data);
void blade_ata_c_register_input_buffer_prefetch_cb(blade_stateful_cb* f);
void blade_ata_c_register_input_buffer_fetch_cb(blade_input_buffer_fetch_cb* f);
void blade_ata_c_register_input_buffer_enqueued_cb(blade_input_buffer_enqueued_cb* f);
void blade_ata_c_register_input_buffer_ready_cb(blade_input_buffer_ready_cb* f);
void blade_ata_c_register_output_buffer_fetch_cb(blade_output_buffer_fetch_cb* f);
void blade_ata_c_register_output_buffer_ready_cb(blade_output_buffer_ready_cb* f);

void blade_ata_c_register_blade_queued_input_clear_cb(blade_clear_queued_cb* f);
void blade_ata_c_register_blade_queued_output_clear_cb(blade_clear_queued_cb* f);

bool blade_ata_c_compute_step();

void blade_ata_c_terminate();

#endif
