#ifndef BLADE_ATA_MODE_h_H
#define BLADE_ATA_MODE_h_H

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <complex.h>

#include "hpguppi_blade_capi.h"
#include "hpguppi_blade_ata_mode_h_config.h"

struct blade_ata_mode_h_config {
    uint32_t channelizerRate;
    uint32_t beamformerBeams;
    uint32_t accumulateRate;
    uint32_t numberOfOutputPolarizations;
    uint32_t integrationSize;

    uint32_t castBlockSize;
    uint32_t channelizerBlockSize;
    uint32_t beamformerBlockSize;
};

static const struct blade_ata_mode_h_config BLADE_ATA_MODE_H_CONFIG = {
    BLADE_ATA_MODE_H_CHANNELIZER_RATE, // .channelizerRate
    BLADE_ATA_MODE_H_OUTPUT_NBEAM, // .beamformerBeams
    BLADE_ATA_MODE_H_ACCUMULATE_RATE, // .acumulateRate
    BLADE_ATA_MODE_H_OUTPUT_NPOL, // .numberOfOutputPolarizations
    BLADE_ATA_MODE_H_INTEGRATION_SIZE, // .integrationSize

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

void blade_ata_h_set_block_time_mjd(double mjd);
void blade_ata_h_set_block_dut1(double dut1);

size_t blade_ata_h_accumulator_counter();

#endif
