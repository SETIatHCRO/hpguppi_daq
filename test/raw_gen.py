from astropy import units as u
import numpy as np
import setigen as stg

n_pols = 2
n_bits = 8
n_ant = 1
n_chan_perant = 1024
n_time = 30720

sample_rate = 2.048e9

file_bytesize = 18*(2**30)
block_bytesize = (n_time*n_ant*n_chan_perant*n_pols*n_bits)//8

n_obschan = n_ant*n_chan_perant
blocks_per_file = file_bytesize // block_bytesize


delays = np.array([int(i*5e-8*sample_rate) for i in range(n_ant)]) # as samples
antenna_array = stg.voltage.MultiAntennaArray(num_antennas=n_ant,
                                    sample_rate=sample_rate,
                                    fch1=6*u.GHz,
                                    ascending=False,
                                    num_pols=2,
                                    delays=delays)

for stream in antenna_array.bg_streams:
        stream.add_noise(v_mean=0,
                        v_std=1)

        stream.add_constant_signal(f_start=6002.2e6*u.Hz,
                                drift_rate=-2*u.Hz/u.s,
                                level=0.002)

digitizer = stg.voltage.RealQuantizer(target_fwhm=32,
                                      num_bits=8)

filterbank = stg.voltage.PolyphaseFilterbank(num_taps=8,
                                             num_branches=n_chan_perant*2)

requantizer = stg.voltage.ComplexQuantizer(target_fwhm=32,
                                           num_bits=n_bits)

rvb = stg.voltage.RawVoltageBackend(antenna_array,
                                    digitizer=digitizer,
                                    filterbank=filterbank,
                                    requantizer=requantizer,
                                    start_chan=0,
                                    num_chans=n_chan_perant,
                                    block_size=block_bytesize,
                                    blocks_per_file=blocks_per_file,
                                    num_subblocks=32)

rvb.record(
	output_file_stem='./synth_large',
	num_blocks=(2**34)//block_bytesize,
	length_mode='num_blocks',
	header_dict={
        'TELESCOP': 'SETIGEN',
        'OBSID'   : 'SYNTHETIC',
	},
	verbose=False,
	load_template=False
)