#ifndef BLADE_ATA_MODE_B_CONFIG_H
#define BLADE_ATA_MODE_B_CONFIG_H

#include "config_input_shape.h"

#define BLADE_ATA_MODE_B_CIRCULAR_POLARIZATION false
#define BLADE_ATA_MODE_B_CHANNELIZER_RATE 1 // [1, 4]; <= 1 mitigates the channlisation

#define BLADE_ATA_MODE_B_OUTPUT_NBEAM 2
#define BLADE_ATA_MODE_B_OUTPUT_NCOMPLEX_BYTES 4

#if BLADE_ATA_MODE_B_OUTPUT_NCOMPLEX_BYTES == 8
	#define BLADE_ATA_MODE_B_OUTPUT_ELEMENT_T CF32
#else
	#define BLADE_ATA_MODE_B_OUTPUT_ELEMENT_T CF16
#endif

#define BLADE_ATA_MODE_B_DATA_SIZE (BLADE_ATA_MODE_B_OUTPUT_NBEAM *\
                               N_INPUT_CHANNELS *\
                               N_INPUT_BLOCK_TIME *\
                               N_INPUT_POL *\
                               BLADE_ATA_MODE_B_OUTPUT_NCOMPLEX_BYTES)

#endif