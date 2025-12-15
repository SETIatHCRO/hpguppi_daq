#ifndef BLADE_ATA_MODE_B_H
#define BLADE_ATA_MODE_B_H

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <complex.h>

#include "hpguppi_blade_capi.h"
#include "hpguppi_blade_ata_mode_b_config.h"

struct blade_ata_mode_b_config {
    uint32_t channelizerRate;
    uint32_t beamformerBeams;

    uint32_t castBlockSize;
    uint32_t channelizerBlockSize;
    uint32_t beamformerBlockSize;
};

static const struct blade_ata_mode_b_config BLADE_ATA_MODE_B_CONFIG = {
    BLADE_ATA_MODE_B_CHANNELIZER_RATE, // .channelizerRate
    BLADE_ATA_MODE_B_OUTPUT_NBEAM, // .beamformerBeams

    512, // .castBlockSize
    512, // .channelizerBlockSize
    512  // .beamformerBlockSize
};

bool blade_ata_b_initialize(
    struct blade_ata_mode_b_config ata_b_config,
    size_t numberOfWorkers,
    struct blade_ata_observation_meta* observationMeta,
    struct LonLatAlt* arrayReferencePosition,
    double* obs_phase_center_radecrad,
    double* beamCoordinates_radecrad,
    double* antennaPositions_xyz,
    double _Complex* antennaCalibrations
);

void blade_ata_b_set_block_time_mjd(double mjd);
void blade_ata_b_set_block_dut1(double dut1);

#endif
