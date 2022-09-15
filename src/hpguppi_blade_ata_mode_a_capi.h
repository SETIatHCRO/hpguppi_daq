#ifndef BLADE_ATA_MODE_A_H
#define BLADE_ATA_MODE_A_H

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <complex.h>

#include "hpguppi_blade_capi.h"
#include "hpguppi_blade_ata_structs.h"
#include "hpguppi_blade_ata_mode_a_config.h"

struct blade_ata_mode_a_config {
    struct blade_ata_input_dims inputDims;
    uint32_t channelizerRate;
    uint32_t beamformerBeams;
    uint32_t integrationSize;
    uint32_t numberOfOutputPolarizations;

    uint32_t outputMemWidth;
    uint32_t outputMemPad;

    uint32_t castBlockSize;
    uint32_t channelizerBlockSize;
    uint32_t beamformerBlockSize;
    uint32_t detectorBlockSize;
};

static const struct blade_ata_mode_a_config BLADE_ATA_MODE_A_CONFIG = {
    { // .inputDims
        1, // .NBEAMS
        BLADE_ATA_MODE_A_INPUT_NANT, // .NANTS 
        BLADE_ATA_MODE_A_ANT_NCHAN, // .NCHANS
        BLADE_ATA_MODE_A_NTIME, // .NTIME 
        BLADE_ATA_MODE_A_NPOL, // .NPOLS 
    }, // .inputDims
    BLADE_ATA_MODE_A_CHANNELIZER_RATE, // .channelizerRate
    BLADE_ATA_MODE_A_OUTPUT_NBEAM, // .beamformerBeams
    BLADE_ATA_MODE_A_INTEGRATION_SIZE, // .integrationSize
    BLADE_ATA_MODE_A_OUTPUT_NPOL, // .numberOfOutputPolarizations

    BLADE_ATA_MODE_A_OUTPUT_MEMCPY2D_WIDTH, // .outputMemWidth
    BLADE_ATA_MODE_A_OUTPUT_MEMCPY2D_PAD, // .outputMemPad

    512, // .castBlockSize
    512, // .channelizerBlockSize
    512, // .beamformerBlockSize
    512  // .detectorBlockSize
};

bool blade_ata_a_initialize(
    struct blade_ata_mode_a_config ata_a_config,
    size_t numberOfWorkers,
    struct blade_ata_observation_meta* observationMeta,
    struct LonLatAlt* arrayReferencePosition,
    double* obs_phase_center_radecrad,
    double* beamCoordinates_radecrad,
    double* antennaPositions_xyz,
    double _Complex* antennaCalibrations
);
size_t blade_ata_a_get_output_size();

void blade_ata_a_set_block_time_mjd(double mjd);
void blade_ata_a_set_block_dut1(double dut1);

void blade_ata_a_register_user_data(void* user_data);
void blade_ata_a_register_input_buffer_prefetch_cb(blade_stateful_cb* f);
void blade_ata_a_register_input_buffer_fetch_cb(blade_input_buffer_fetch_cb* f);
void blade_ata_a_register_input_buffer_enqueued_cb(blade_input_buffer_enqueued_cb* f);
void blade_ata_a_register_input_buffer_ready_cb(blade_input_buffer_ready_cb* f);
void blade_ata_a_register_output_buffer_fetch_cb(blade_output_buffer_fetch_cb* f);
void blade_ata_a_register_output_buffer_ready_cb(blade_output_buffer_ready_cb* f);

void blade_ata_a_register_blade_queued_input_clear_cb(blade_clear_queued_cb* f);
void blade_ata_a_register_blade_queued_output_clear_cb(blade_clear_queued_cb* f);

bool blade_ata_a_compute_step();

void blade_ata_a_terminate();

#endif
