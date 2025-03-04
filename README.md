# Hashpipe SETI Data Acquisition

----------------------------------------------

# Dockerfile usage

`docker build -t hp_daq .`

`docker run -it --cap-add IPC_LOCK --volume ${PWD}/test:/mnt/buf0 --env HASHPIPE_KEYFILE=/work --env LD_LIBRARY_PATH=/work/hashpipe/src/.libs:/work/xGPU/src:/work/uvh5c99/build:/work/uvh5c99/build/subprojects/radiointerferometryc99:${LD_LIBRARY_PATH} hp_daq`

# Basic library dependencies

## [Hashpipe](https://github.com/david-macmahon/hashpipe/)

This is foundational of course.

```
$ git clone -b seti https://github.com/MydonSolutions/hashpipe
$ cd hashpipe
$ git checkout 81a79e626d4fe78f3f7cc6209be45b8569fae42d
$ autoreconf -is && ./configure && make
```

## [Filterbankh5c99](https://github.com/MydonSolutions/filterbankh5c99)

Hpguppi_daq currently expects this for certain files that aren't always in use:
- hpguppi_ata_fildisk_thread # Contemporary filterbank output

```
$ git clone https://github.com/MydonSolutions/filterbankh5c99
$ cd filterbankh5
$ git submodule init
$ meson build -Dprefix=${PWD}/install
$ cd build
$ ninja test && ninja install
```

[CUDA](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html) is required.

## [SLA](https://github.com/scottransom/pyslalib)
- [ ] ? Absorb into radiointerferometryc99.

Hpguppi_daq uses in only a few critical places:
- hpguppi_params # manage struct meta-data from status key-values
- hpguppi_time # MJD time calculation

```
$ git clone https://github.com/scottransom/pyslalib
$ cd pyslalib
$ make libsla.so
```

# Situational library dependencies

## [UVH5C99](https://github.com/MydonSolutions/uvh5c99)

This library simplifies the authoring of UVH5 files in C. It's a dependency of:
- hpguppi_ata_obs_uvh5disk_thread

```
$ git clone https://github.com/MydonSolutions/uvh5c99
$ cd uvh5c99
$ git submodule init
$ meson build
$ cd build
$ ninja test
```

## [XGPU](https://github.com/GPU-correlators/xGPU)

Correlation is accomplished using xGPU. It's a dependency of:
- hpguppi_ata_xgpu_thread
- hpguppi_xgpu_databuf
- hpguppi_atasnap_xgpu_disk_thread

```
$ git clone https://github.com/GPU-correlators/xGPU
$ cd xGPU/src
$ make clean && make NTIME=32768 NTIME_PIPE=128 NPOL=2 NFREQUENCY=512 NSTATION=16 CUDA_ARCH=sm_86 DP4A=yes
```

## [BLADE](https://github.com/luigifcruz/blade)

Beamforming is accomplished with BLADE (^0.7.0) empowered threads:
- hpguppi_ata_blade_beamformer_thread

```
$ git clone -b blade-cli https://github.com/luigifcruz/blade
$ cd blade
$ CC=gcc-10 CXX=g++-10 meson build -Dprefix=${PWD}/install
$ cd build
$ ninja test && ninja install
```

# Compilation

The `$ ./configure` step determines which threads are compiled and indeed available for deployment.

- [`--with-sla-lib`](#slahttpsgithubcomscottransompyslalib) is required for MJD time calculations. 
Provide the directory containing the library.

- [`--with-blade`](#bladehttpsgithubcomluigifcruzblade) will enable BLADE related threads. 
For this `CXX=g++-10` should be set.
Provide the installation directory of BLADE (-Dprefix defaults to /usr/local).

- `--with-cuda-include` required to compile BLADE related C-API.

- [`--with-xgpu`](https://github.com/GPU-correlators/xGPU) will enable xGPU related threads.
Provide the compilation directory (i.e. xGPU/src).

- [`--with-uvh5`](https://github.com/MydonSolutions/uvh5c99) will enable UVH5 related threads.
Provide the compilation directory (i.e. uvh5c99/build).

- [`--with-filterbankh5c99`](https://github.com/MydonSolutions/filterbankh5c99) will enable Filterbank related threads.

All in all:
`hpguppi_daq/src$ CXX=g++-11 ./configure --with-libsla=../../pyslalib --with-hashpipe=../../hashpipe/src/.libs --with-cuda-include=/usr/local/cuda-11.4.1/include --with-xgpu=../../xGPU/src --with-uvh5=../../uvh5c99/build --with-rawspec=../../rawspec --with-blade=../../blade/install --with-filterbankh5c99=../../filterbankh5c99/install

Thereafter, `$ make` compiles the threads into `$ ./.libs/hpguppi_daq.so`.

# Instantiation

The `init_hpguppi.py` script assesses the `config_hpguppi.yaml` file to achieve the instance-system named as the first positional argument. The script is typically executed with super-user privileges.

# Further Notes

## [rb-hashpipe](https://github.com/david-macmahon/rb-hashpipe)
```
$ git clone https://github.com/david-macmahon/rb-hashpipe
$ cd rb-hashpipe
$ rake package 
$ sudo apt install libncurses5-dev
$ cd pkg
$ sudo gem install --local ./hashpipe-0.6.3.gem -- --with-hashpipe-include=../../hashpipe/src --with-hashpipestatus-lib=../../hashpipe/src/.libs
```

## System optimisation

[Optimal BIOS settings](https://hpcadvisorycouncil.atlassian.net/wiki/spaces/HPCWORKS/pages/1280442391/AMD+2nd+Gen+EPYC+CPU+Tuning+Guide+for+InfiniBand+HPC?focusedCommentId=2152333319)
- APBDIS = 1 [NB Configuration]
- Fixed SOC Pstate = P0 [NB Configuration]
- IOMMU = disabled [NB Configuration]
- DF Cstates = disabled [NB Configuration]
- SMT = disabled [CPU Configuration]
- Local APIC mode = x2APIC [CPU Configuration]
- NUMA nodes per socket = NPS1 [ACPI Settings]
- L3 cache as NUMA = disabled [ACPI Settings (ACPI SRAT L3 Cache As NUMA Domain)] {{{ it’s recommended to be enabled, but then causes too many NUMA nodes }}}
- PCIe Relaxed Ordering = enabled [PCIe/PCI/PnP Configuration]
