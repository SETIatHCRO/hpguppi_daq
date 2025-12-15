#ifndef BLADE_ATA_MODE_X_H
#define BLADE_ATA_MODE_X_H

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <complex.h>

#include "hpguppi_blade_capi.h"
#include "hpguppi_blade_ata_mode_x_config.h"

struct blade_ata_mode_x_config {
    uint32_t channelizerRate;
    uint32_t beamformerBeams;
    uint32_t integrationSize;
    uint32_t frequencyIntegrationSize;
    bool kurtosisEnabled;
    int32_t kurtosisSigma;
    char kurtosisMaskOutputFilepath[256];

    uint32_t castBlockSize;
    uint32_t channelizerBlockSize;
    uint32_t correlatorBlockSize;
};

static const struct blade_ata_mode_x_config BLADE_ATA_MODE_X_CONFIG = {
    BLADE_ATA_MODE_X_CHANNELIZER_RATE, // .channelizerRate
    1, // .beamformerBeams
    BLADE_ATA_MODE_X_INTEGRATION_SIZE, // .integrationSize
    1, // .frequencyIntegrationSize

    false, // .kurtosisEnabled
    3, // .kurtosisSigma
    "/dev/null", // .kurtosisMaskOutputFilepath

    512, // .castBlockSize
    512, // .channelizerBlockSize
    32, // .correlatorBlockSize
};

bool blade_ata_x_initialize(
    struct blade_ata_mode_x_config ata_x_config
);

#endif
