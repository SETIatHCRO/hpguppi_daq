#ifndef BLADE_ATA_MODE_X_CONFIG_H
#define BLADE_ATA_MODE_X_CONFIG_H

#include "config_input_shape.h"

#define BLADE_ATA_MODE_X_INTEGRATION_SIZE 262144 // just a fallback default
#define BLADE_ATA_MODE_X_INTEGRATION_FACTOR 1
#define BLADE_ATA_MODE_X_CHANNELIZER_RATE 1 // 1 mitigates the channelization

#if BLADE_ATA_MODE_X_CHANNELIZER_RATE == 1
// Continuum correlator
#else
// spectral-line correlator # 16 MHz
#endif

#define BLADE_ATA_MODE_X_CONJUGATION_INDEX 0

#define BLADE_ATA_MODE_X_OUTPUT_ELEMENT_T CF32
#define BLADE_ATA_MODE_X_OUTPUT_NCOMPLEX_BYTES 8

#define BLADE_ATA_MODE_X_DATA_SIZE ((size_t)(N_INPUT_ASPECTS*(N_INPUT_ASPECTS+1)/2) *\
                               (size_t)(N_INPUT_CHANNELS*BLADE_ATA_MODE_X_CHANNELIZER_RATE) *\
                               1 *\
                               N_INPUT_POL*N_INPUT_POL *\
                               BLADE_ATA_MODE_X_OUTPUT_NCOMPLEX_BYTES)

#endif