####v Nvidia builder v####
FROM nvidia/cuda:12.2.0-devel-ubuntu22.04 AS nvidia_builder
# RUN apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/3bf863cc.pub
# RUN  apt-get update --fix-missing

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get --fix-missing update -y && apt-get install -y \
    g++-10 \
    gcc-10 \
    git \
    liberfa-dev \
    libfmt-dev \
    libhdf5-dev \
    libboost-all-dev \
    libbenchmark-dev \
    libgtest-dev \
    libspdlog-dev \
    build-essential \
    # GENERAL
    numactl \
    linux-tools-generic \
    # HPGUPPI_DAQ
    libcfitsio-dev \
    # PYSLALIB
    gfortran \
    # RB-HASHPIPE
    ruby-dev \
    libncurses5-dev \
    pkg-config \
    cmake \
    python3-dev \
    python3-pip
    # libpq-dev

RUN python3 -m pip install meson ninja numpy astropy pandas
####^ Nvidia builder ^####

####v BLADE Builder v####
FROM nvidia_builder AS blade_builder

ENV CUDA_ROOT=/usr/local/cuda
WORKDIR /work

RUN cd /work \
&& git clone https://github.com/luigifcruz/blade \
&& cd blade \
&& git submodule update --init --recursive \
&& CC=gcc-10 CXX=g++-10 meson setup build -Dprefix=${PWD}/install \
&& cd build \
&& ninja install
####^ BLADE Builder ^####

####v Hashpipe Builder v####
FROM nvidia_builder AS hashpipe_builder

WORKDIR /work

RUN cd /work \
&& git clone -b seti https://github.com/MydonSolutions/hashpipe \
&& cd hashpipe/src \
&& git checkout 81a79e626d4fe78f3f7cc6209be45b8569fae42d \
&& autoreconf -is \
&& ./configure \
&& make
####^ Hashpipe Builder ^####

####v SLA Builder v####
FROM nvidia_builder AS sla_builder

WORKDIR /work

RUN cd /work \
&& git clone https://github.com/scottransom/pyslalib &&\
    cd pyslalib \
&& make libsla.so
####^ SLA Builder ^####

# ####v xGPU Builder v####
# FROM nvidia_builder AS xgpu_builder

# WORKDIR /work

# RUN cd /work \
# && git clone https://github.com/GPU-correlators/xGPU \
# && cd xGPU/src \
# && make clean \
# && make NTIME=32768 NTIME_PIPE=128 NPOL=2 NFREQUENCY=512 NSTATION=16 CUDA_ARCH=sm_86 DP4A=yes
# ####^ xGPU Builder ^####

####v UVH5C99 Builder v####
FROM nvidia_builder AS uvh5c99_builder

WORKDIR /work

RUN cd /work \
&& git clone https://github.com/MydonSolutions/uvh5c99 \
&& cd uvh5c99 \
&& git submodule update --init \
&& meson setup build -Dprefix=${PWD}/install \
&& cd build \
&& ninja install
####^ UVH5C99 Builder ^####

####v FILTERBANKC99 Builder v####
FROM nvidia_builder AS filterbankc99_builder

WORKDIR /work

RUN cd /work \
&& git clone https://github.com/MydonSolutions/filterbankc99 \
&& cd filterbankc99 \
&& git submodule update --init \
&& meson setup build -Dprefix=${PWD}/install \
&& cd build \
&& ninja install
####^ FILTERBANKC99 Builder ^####

####v HPDAQ Builder v####
FROM nvidia_builder as hpdaq_builder

# yaml for the init_hpguppi.py script
# RUN python3 -m pip install yamlpy

# Install dependencies
WORKDIR /work
RUN mkdir /work/logs

## Hashpipe
COPY --from=hashpipe_builder /work/hashpipe /work/hashpipe

# ## RB-HASHPIPE
# RUN cd /work \
# && gem install redis \
# && git clone https://github.com/david-macmahon/rb-hashpipe \
# && cd rb-hashpipe \
# && rake package \
# && cd pkg \
# && gem install \
#     --local ./hashpipe-0.6.3.gem -- \
#     --with-hashpipe-include=/work/hashpipe/src \
#     --with-hashpipestatus-lib=/work/hashpipe/src/.libs \
# && gem install curses

## BLADE
COPY --from=blade_builder /work/blade /work/blade

## SLA
COPY --from=sla_builder /work/pyslalib /work/pyslalib

## xGPU
# COPY --from=xgpu_builder /work/xGPU /work/xGPU

## UVH5C99
COPY --from=uvh5c99_builder /work/uvh5c99 /work/uvh5c99

## FILTERBANKC99
COPY --from=filterbankc99_builder /work/filterbankc99 /work/filterbankc99

## Hpguppi_daq
COPY . /work/hpguppi_daq
RUN cd /work/hpguppi_daq/src \
&& git submodule update --init \
&& cd ata_antenna_weights_binary && meson setup build && cd build && ninja install && cd ../../ \
&& rm -rf install-sh libtool ltmain.sh install-sh depcomp configure config.* aclocal.m4 compile automa4te.cache/ missing Makefile Makefile.in \
&& autoreconf -is \
&& mkdir ../build \
&& cd ../build \
&& CXX=g++-10 ../src/configure \
    --with-sla-lib=/work/pyslalib \
    --with-hashpipe=/work/hashpipe/src/.libs \
    --with-cuda-include=/usr/local/cuda-12.2/include \
    --with-filterbankc99=/work/filterbankc99/install \
    --with-uvh5=/work/uvh5c99/install \
    --with-blade=/work/blade/install
# && make
    # --with-xgpu=/work/xGPU/src 
####^ HPDAQ Builder ^####