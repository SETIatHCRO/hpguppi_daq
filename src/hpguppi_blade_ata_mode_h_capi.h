#ifndef BLADE_ATA_MODE_h_H
#define BLADE_ATA_MODE_h_H

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <complex.h>

#include "hpguppi_blade_capi.h"
#include "hpguppi_blade_ata_structs.h"
#include "hpguppi_blade_ata_mode_h_config.h"

struct blade_ata_mode_h_config {
    struct blade_ata_input_dims inputDims;
    uint32_t channelizerRate;
    uint32_t beamformerBeams;
    uint32_t accumulateRate;
    uint32_t numberOfOutputPolarizations;
    uint32_t integrationSize;

    uint32_t outputMemWidth;
    uint32_t outputMemPad;

    uint32_t castBlockSize;
    uint32_t channelizerBlockSize;
    uint32_t beamformerBlockSize;
};

static const struct blade_ata_mode_h_config BLADE_ATA_MODE_H_CONFIG = {
    { // .inputDims
        1, // .NBEAMS
        BLADE_ATA_MODE_H_INPUT_NANT, // .NANTS 
        BLADE_ATA_MODE_H_ANT_NCHAN, // .NCHANS
        BLADE_ATA_MODE_H_NTIME, // .NTIME 
        BLADE_ATA_MODE_H_NPOL, // .NPOLS 
    }, // .inputDims
    BLADE_ATA_MODE_H_CHANNELIZER_RATE, // .channelizerRate
    BLADE_ATA_MODE_H_OUTPUT_NBEAM, // .beamformerBeams
    BLADE_ATA_MODE_H_ACCUMULATE_RATE, // .acumulateRate
    BLADE_ATA_MODE_H_OUTPUT_NPOL, // .numberOfOutputPolarizations
    1, // .integrationSize

    BLADE_ATA_MODE_H_OUTPUT_MEMCPY2D_WIDTH, // .outputMemWidth
    BLADE_ATA_MODE_H_OUTPUT_MEMCPY2D_PAD, // .outputMemPad

    512, // .castBlockSize
    512, // .channelizerBlockSize
    512  // .beamformerBlockSize
};

bool blade_ata_h_initialize(
    struct blade_ata_mode_h_config ata_h_config,
    size_t numberOfWorkers,
    struct blade_ata_observation_meta* observationMeta,
    struct LonLatAlt* arrayReferencePosition,
    double* obs_phase_center_radecrad,
    double* beamCoordinates_radecrad,
    double* antennaPositions_xyz,
    double _Complex* antennaCalibrations
);
size_t blade_ata_h_get_input_size();
size_t blade_ata_h_get_output_size();

size_t blade_ata_h_accumulator_counter();

bool blade_ata_h_enqueue_b(void* input_ptr, const size_t b_id, double time_mjd, double dut1);
bool blade_ata_h_dequeue_b(size_t* b_id);
bool blade_ata_h_enqueue_h(void* output_ptr, const size_t h_id);
bool blade_ata_h_dequeue_h(size_t* id);

void blade_ata_h_terminate();

#endif
