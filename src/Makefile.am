ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS     =

AM_CPPFLAGS += -I"@HASHPIPE_INCDIR@"

# AM_CFLAGS is used for all C compiles
AM_CFLAGS = -fopenmp -ggdb -fPIC -O3 -Wall -Werror -Wno-error=cpp -fno-strict-aliasing @X86_FEATURE_CFLAGS@
AM_CFLAGS += -Wno-address-of-packed-member 
AM_CFLAGS += -Wno-unused-const-variable
AM_CFLAGS += -Wno-stringop-truncation
AM_CFLAGS += -Wno-unused-result
AM_CFLAGS += -Wno-format-truncation
AM_CFLAGS += -Wno-unused-but-set-variable

# Convenience variables to group source files
hpguppi_databuf = hpguppi_databuf.h \
		  hpguppi_databuf.c

hpguppi_support = hpguppi_atasnap.h \
		  hpguppi_atasnap.c \
		  hpguppi_params.c \
		  hpguppi_pktbuf.h \
		  hpguppi_pktbuf.c \
		  hpguppi_time.h   \
		  hpguppi_time.c   \
		  hpguppi_util.h   \
		  hpguppi_util.c   \
		  polyco_struct.h  \
		  psrfits.h

# hpguppi_support += \
#		  hpguppi_udp.h    \
#		  hpguppi_udp.c

libsla_support = slalib.h sla.c f77.h

hpguppi_threads = \
		hpguppi_ata_ibv_payload_order_thread.c \
		hpguppi_atasnap_obs_rawdisk_thread.c
#		guppiraw_ibv_packetize_thread.c # depends on rawspec currently

if IBVERBS_ENABLED
hpguppi_daq_la_LIBADD = -lhashpipe_ibverbs
hpguppi_threads += hpguppi_ibverbs_pkt_thread.c
else
hpguppi_daq_la_LIBADD = -lhashpipe
endif

# This is the hpguppi_daq plugin
lib_LTLIBRARIES = hpguppi_daq.la
hpguppi_daq_la_SOURCES  = $(hpguppi_databuf) $(hpguppi_support) $(libsla_support) $(hpguppi_threads)
hpguppi_daq_la_LIBADD  += -lsla -lrt
hpguppi_daq_la_CFLAGS 	= $(AM_CFLAGS)
hpguppi_daq_la_LDFLAGS  = -avoid-version -module -shared -export-dynamic
hpguppi_daq_la_LDFLAGS += -L"@LIBSLADIR@" -Wl,-rpath,"@LIBSLADIR@"
hpguppi_daq_la_LDFLAGS += -L"@HASHPIPE_LIBDIR@" -Wl,-rpath,"@HASHPIPE_LIBDIR@"

if BLADE_ENABLED
# bladeapi library
AM_CPPFLAGS += -I"@BLADE_INCDIR@"
AM_CFLAGS += -I"ata_antenna_weights_binary/include"

hpguppi_support += \
		ata_antenna_weights_binary/src/antenna_weights.c \
		ata_antenna_weights_binary/include/antenna_weights.h

hpguppi_blade_c_api_files = \
		hpguppi_blade_capi.h \
		hpguppi_blade_capi.cpp \
		hpguppi_blade_ata_structs.h \
		hpguppi_blade_ata_mode.h \
		hpguppi_blade_ata_mode_a_capi.cpp \
		hpguppi_blade_ata_mode_a_capi.h \
		hpguppi_blade_ata_mode_a_config.h \
		hpguppi_blade_ata_mode_h_capi.cpp \
		hpguppi_blade_ata_mode_h_capi.h \
		hpguppi_blade_ata_mode_h_config.h \
		hpguppi_blade_ata_mode_b_capi.cpp \
		hpguppi_blade_ata_mode_b_capi.h \
		hpguppi_blade_ata_mode_b_config.h

noinst_LTLIBRARIES = libbladeapi.la
libbladeapi_la_SOURCES = $(hpguppi_blade_c_api_files)
libbladeapi_la_CXXFLAGS = -std=c++2a -I "@BLADE_INCDIR@"
libbladeapi_la_LDFLAGS = -L"@BLADE_LIBDIR@"
libbladeapi_la_LIBADD = -lblade
if CUDA_ENABLED
libbladeapi_la_CXXFLAGS  += -I "@CUDA_INCDIR@"
endif

# hpguppi extension
hpguppi_support += \
		hpguppi_blade_ata_mode_b_capi.h \
		hpguppi_blade_ata_mode_b_config.h
hpguppi_databuf += \
		hpguppi_blade_databuf.c \
		hpguppi_blade_databuf.h

hpguppi_threads += \
		hpguppi_ata_blade_beamformer_thread.c

hpguppi_daq_la_LIBADD  += -lblade libbladeapi.la
hpguppi_daq_la_LDFLAGS += -L"@BLADE_LIBDIR@" -Wl,-rpath,"@BLADE_LIBDIR@"
hpguppi_daq_la_CFLAGS  += -I "@BLADE_INCDIR@"
endif

if FILTERBANKC99_ENABLED
hpguppi_daq_la_LDFLAGS += -L"@FILTERBANKC99_LIBDIR@" -Wl,-rpath,"@FILTERBANKC99_LIBDIR@"
hpguppi_daq_la_LIBADD  += -lfilterbankc99
hpguppi_daq_la_LIBADD  += -lh5dsc99
hpguppi_daq_la_CFLAGS  += -I "@FILTERBANKC99_INCDIR@"

if BLADE_ENABLED
hpguppi_threads += \
		hpguppi_ata_fildisk_thread.c
endif
endif

if CUDA_ENABLED
hpguppi_daq_la_CFLAGS  += -I "@CUDA_INCDIR@"
endif

if XGPU_ENABLED
hpguppi_databuf += \
		hpguppi_xgpu_databuf.c \
		hpguppi_xgpu_databuf.h

hpguppi_threads += \
		hpguppi_ata_xgpu_thread.c
if UVH5C99_ENABLED
hpguppi_threads += \
		hpguppi_ata_xintegration_thread.c
endif

hpguppi_daq_la_LDFLAGS += -L"@XGPU_LIBDIR@" -lxgpu -Xlinker -lrt
hpguppi_daq_la_CFLAGS  += -I "@XGPU_INCDIR@"

if UVH5C99_ENABLED
hpguppi_threads += \
		hpguppi_ata_obs_uvh5disk_thread.c
hpguppi_daq_la_LDFLAGS += -L"@UVH5C99_LIBDIR@"
hpguppi_daq_la_LIBADD += -luvh5
hpguppi_daq_la_CFLAGS  += -I "@UVH5C99_INCDIR@"
hpguppi_daq_la_CFLAGS  += -I "@UVH5C99_INCDIR@/../external_src/tomlc99"
endif # UVH5C99_ENABLED 

endif


if RADIOINTERFEROMETRYC99_ENABLED
hpguppi_daq_la_LDFLAGS += -L"@RADIOINTERFEROMETRYC99_LIBDIR@"
hpguppi_daq_la_LIBADD += -lradiointerferometryc99
hpguppi_daq_la_CFLAGS  += -I "@RADIOINTERFEROMETRYC99_INCDIR@"
endif # RADIOINTERFEROMETRYC99_ENABLED 

# Installed scripts
dist_bin_SCRIPTS = init_hpguppi.py


# vi: set ts=8 noet :
