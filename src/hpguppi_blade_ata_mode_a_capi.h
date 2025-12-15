#ifndef BLADE_ATA_MODE_A_H
#define BLADE_ATA_MODE_A_H

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <complex.h>

#include "hpguppi_blade_capi.h"
#include "hpguppi_blade_ata_mode_a_config.h"

struct blade_ata_mode_a_config {
    uint32_t channelizerRate;
    uint32_t beamformerBeams;
    uint32_t integrationSize;
    uint32_t numberOfOutputPolarizations;

    bool kurtosisEnabled;
    uint32_t kurtosisSigma;
    uint32_t kurtosisChannelLength;
    uint32_t kurtosisNumberOfMaskRuns;
    char kurtosisMaskOutputFilepath[256];

    uint32_t castBlockSize;
    uint32_t channelizerBlockSize;
    uint32_t beamformerBlockSize;
    uint32_t detectorBlockSize;
};

static const struct blade_ata_mode_a_config BLADE_ATA_MODE_A_CONFIG = {
    BLADE_ATA_MODE_A_CHANNELIZER_RATE, // .channelizerRate
    BLADE_ATA_MODE_A_OUTPUT_NBEAM, // .beamformerBeams
    BLADE_ATA_MODE_A_INTEGRATION_SIZE, // .integrationSize
    BLADE_ATA_MODE_A_OUTPUT_NPOL, // .numberOfOutputPolarizations

    false, // .kurtosisEnabled
    3, // .kurtosisSigma
    256, // .kurtosisChannelLength
    128, // .kurtosisNumberOfMaskRuns
    "/dev/null", // .kurtosisMaskOutputFilepath

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

void blade_ata_a_set_block_time_mjd(double mjd);
void blade_ata_a_set_block_dut1(double dut1);

#endif
