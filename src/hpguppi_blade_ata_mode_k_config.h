#ifndef BLADE_ATA_MODE_K_CONFIG_H
#define BLADE_ATA_MODE_K_CONFIG_H

#define BLADE_ATA_MODE_K_OUTPUT_ELEMENT_T CF32
#define BLADE_ATA_MODE_K_OUTPUT_NCOMPLEX_BYTES 8

#define BLADE_ATA_MODE_K_DATA_SIZE ((size_t)N_INPUT_ASPECTS *\
                               N_INPUT_CHANNELS *\
                               N_INPUT_BLOCK_TIME *\
                               N_INPUT_POL *\
                               BLADE_ATA_MODE_K_OUTPUT_NCOMPLEX_BYTES)

#endif