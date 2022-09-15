#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <math.h>

#include <complex.h>

#include "hashpipe.h"
#include "hpguppi_time.h"
#include "hpguppi_databuf.h"
#include "hpguppi_ata_blade_mode.h"
#include "hpguppi_blade_databuf.h"

#include "uvh5/uvh5_toml.h"
#include "radiointerferometryc99.h"
#include "antenna_weights.h"

typedef struct {
    const int blade_number_of_workers;

    hashpipe_status_t* status;

    const char* thread_name;
    const char* status_key;
    char status_state;
    struct timespec ts_last_status_update;
    struct timespec ts_buffer_wait_timeout;

    int32_t prev_flagged_NANTS;
    int32_t prev_flagged_NCHAN;
    int32_t prev_flagged_NTIME;
    int32_t prev_flagged_NPOLS;
    int64_t prev_pktidx_obs_start;

    size_t accumulator_counter;

    size_t in_index;
    size_t out_index_free;
    size_t out_index_fill;

    hpguppi_input_databuf_t* in;
    hpguppi_blade_output_databuf_t* out;

    void* out_intermediary[N_INPUT_BLOCKS];

    uint64_t fill_to_free_moving_sum_ns;
    uint64_t fill_to_free_block_ns[N_INPUT_BLOCKS];
    struct timespec ts_blocks_recvd[N_INPUT_BLOCKS];
} blade_userdata_t;

double jd_mid_block(char* databuf_header) {
  uint64_t pktidx, piperblk;
  struct mjd_t mjd = {0};
  uint64_t synctime = 0;
  double chan_bw = 1.0;

  double realtime_secs = 0.0;
  struct timespec ts;


  hgetu8(databuf_header, "PKTIDX", &pktidx);
  hgetu8(databuf_header, "PIPERBLK", &piperblk);
  hgetr8(databuf_header, "CHAN_BW", &chan_bw);
  hgetu8(databuf_header, "SYNCTIME", &synctime);

  // Calc real-time seconds since SYNCTIME for pktidx, taken to be a multiple of PKTNTIME:
  //
  //                          pktidx
  //     realtime_secs = -------------------
  //                        1e6 * chan_bw
  if (chan_bw != 0.0) {
    realtime_secs = (pktidx + piperblk/2) / (1e6 * fabs(chan_bw));
  }

  ts.tv_sec = (time_t)(synctime + rint(realtime_secs));
  ts.tv_nsec = (long)((realtime_secs - rint(realtime_secs)) * 1e9);

  get_mjd_from_timespec(&ts, &(mjd.stt_imjd), &(mjd.stt_smjd), &(mjd.stt_offs));
  //double t = ((double)mjd.stt_imjd) + ((double) mjd.stt_smjd + mjd.stt_offs)/(86400.0);
  //double t2 = ((double)synctime + (double)realtime_secs) / (double) 86400.0 + (double) 40587.0;
  //fprintf(stderr, "diff is: %.16f %.16f\n", t, t2);
  return (((double)mjd.stt_imjd) + ((double) mjd.stt_smjd + mjd.stt_offs)/(86400.0)) + 2400000.5;
}

// Populates `beam_coordinates` which should be nbeams*2 long (RA/DEC_OFFX pairs)
void collect_beamCoordinates(
  int nbeams,
  double* beam_coordinates,
  double* phase_center,
  char* databuf_header
) {

  // Getting phase center
  hgetr8(databuf_header, "RA_STR",  phase_center+0);
  hgetr8(databuf_header, "DEC_STR", phase_center+1);
  phase_center[0] = calc_rad_from_degree(phase_center[0]* 360.0 / 24.0); // convert from hours to degrees
  phase_center[1] = calc_rad_from_degree(phase_center[1]);


  // Getting beam coordinates
  char coordkey[9] = "DEC_OFFX";
  for(int beam_idx = 0; beam_idx < nbeams; beam_idx++) {
    sprintf(coordkey, "RA_OFF%d", beam_idx%10);
    hgetr8(databuf_header, coordkey, beam_coordinates+beam_idx*2+0);
    beam_coordinates[beam_idx*2+0] = calc_rad_from_degree(beam_coordinates[beam_idx*2+0] * 360.0 / 24.0);

    sprintf(coordkey, "DEC_OFF%d", beam_idx%10);
    hgetr8(databuf_header, coordkey, beam_coordinates+beam_idx*2+1);
    beam_coordinates[beam_idx*2+1] = calc_rad_from_degree(beam_coordinates[beam_idx*2+1]);
  }
}

bool blade_cb_input_buffer_prefetch(void* user_data_void) {
  blade_userdata_t* user_data = (blade_userdata_t*) user_data_void;

  struct timespec ts_now;
  char buf_status[80];

  {// poll once for an incoming block
    int hpguppi_databuf_wait_rv = hpguppi_databuf_wait_filled_timeout(
      user_data->in, user_data->in_index,
      &user_data->ts_buffer_wait_timeout
    );

    clock_gettime(CLOCK_MONOTONIC, &ts_now);
    // We perform some status buffer updates every second
    if (ELAPSED_NS(user_data->ts_last_status_update, ts_now) > 1e9) {
      sprintf(buf_status, "%d/%d", hpguppi_databuf_total_status(user_data->out), user_data->in->header.n_block);
      memcpy(&user_data->ts_last_status_update, &ts_now, sizeof(struct timespec));

      hashpipe_status_lock_safe(user_data->status);
      {
        if (user_data->status_state != 1) {
          hputs(user_data->status->buf, user_data->status_key, "inblocked");
          user_data->status_state=1;
        }

        hputi4(user_data->status->buf, "BLDBLKSZ", BLADE_BLOCK_DATA_SIZE);
        hputs(user_data->status->buf, "BLDBUFST", buf_status);
        hputr4(user_data->status->buf, "BLDBLKMS",
          round(
            (double)(user_data->fill_to_free_moving_sum_ns)
            / (user_data->blade_number_of_workers * N_INPUT_BLOCKS)
          ) / 1e6
        );
      }
      hashpipe_status_unlock_safe(user_data->status);
    }

    if (hpguppi_databuf_wait_rv == HASHPIPE_TIMEOUT) {
      return false;
    }
    else if (hpguppi_databuf_wait_rv != HASHPIPE_TIMEOUT && hpguppi_databuf_wait_rv != HASHPIPE_OK) {
      hashpipe_error(user_data->thread_name, "error status_state for input buffer, rv: %i", hpguppi_databuf_wait_rv);
      pthread_exit(NULL);
      return false;
    }
  }

  char* databuf_header = hpguppi_databuf_header(user_data->in, user_data->in_index);

  {// validate apparent dimensions of data
    int32_t input_buffer_dim_NANTS = 0;
    int32_t input_buffer_dim_NCHAN = 0;
    int32_t input_buffer_dim_NTIME = 0;
    int32_t input_buffer_dim_NPOLS = 0;

    char indb_data_dims_good_flag = 1; // innocent until proven guilty
    hgeti4(databuf_header, "NANTS", &input_buffer_dim_NANTS);
    hgeti4(databuf_header, "NCHAN", &input_buffer_dim_NCHAN);
    hgeti4(databuf_header, "PIPERBLK", &input_buffer_dim_NTIME);
    hgeti4(databuf_header, "NPOL", &input_buffer_dim_NPOLS);

    if (input_buffer_dim_NANTS != BLADE_ATA_CONFIG.inputDims.NANTS) {
      indb_data_dims_good_flag = 0;
      if (user_data->prev_flagged_NANTS != input_buffer_dim_NANTS) {
        user_data->prev_flagged_NANTS = input_buffer_dim_NANTS;
        hashpipe_error(user_data->thread_name, "Incoming data_buffer has NANTS %lu != %lu. Ignored.", input_buffer_dim_NANTS, BLADE_ATA_CONFIG.inputDims.NANTS);
      }
    }
    else if (input_buffer_dim_NCHAN != BLADE_ATA_CONFIG.inputDims.NCHANS) {
      indb_data_dims_good_flag = 0;
      if (user_data->prev_flagged_NCHAN != input_buffer_dim_NCHAN) {
        user_data->prev_flagged_NCHAN = input_buffer_dim_NCHAN;
        hashpipe_error(user_data->thread_name, "Incoming data_buffer has NCHANS %lu != %lu. Ignored.", input_buffer_dim_NCHAN, BLADE_ATA_CONFIG.inputDims.NCHANS);
      }
    }
    else if (input_buffer_dim_NTIME != BLADE_ATA_CONFIG.inputDims.NTIME) {
      indb_data_dims_good_flag = 0;
      if (user_data->prev_flagged_NTIME != input_buffer_dim_NTIME) {
        user_data->prev_flagged_NTIME = input_buffer_dim_NTIME;
        hashpipe_error(user_data->thread_name, "Incoming data_buffer has NTIME %lu != %lu. Ignored.", input_buffer_dim_NTIME, BLADE_ATA_CONFIG.inputDims.NTIME);
      }
    }
    else if (input_buffer_dim_NPOLS != BLADE_ATA_CONFIG.inputDims.NPOLS) {
      indb_data_dims_good_flag = 0;
      if (user_data->prev_flagged_NPOLS != input_buffer_dim_NPOLS) {
        user_data->prev_flagged_NPOLS = input_buffer_dim_NPOLS;
        hashpipe_error(user_data->thread_name, "Incoming data_buffer has NPOLS %lu != %lu. Ignored.", input_buffer_dim_NPOLS, BLADE_ATA_CONFIG.inputDims.NPOLS);
      }
    }

    if (!indb_data_dims_good_flag) {
      hpguppi_databuf_set_free(user_data->in, user_data->in_index);
      user_data->in_index  = (user_data->in_index + 1) % user_data->in->header.n_block;
      return false;
    }
    else {
      user_data->prev_flagged_NANTS = input_buffer_dim_NANTS;
      user_data->prev_flagged_NCHAN = input_buffer_dim_NCHAN;
      user_data->prev_flagged_NTIME = input_buffer_dim_NTIME;
      user_data->prev_flagged_NPOLS = input_buffer_dim_NPOLS;
    }
  }

  {// re-setup if a new observation
    int64_t pktidx_obs_start, pktidx_blk_start;
    hgeti8(databuf_header, "PKTSTART", &pktidx_obs_start);
    hgeti8(databuf_header, "BLKSTART", &pktidx_blk_start);

    // if first block of observation
    if (pktidx_obs_start != user_data->prev_pktidx_obs_start && pktidx_obs_start >= pktidx_blk_start) {
      UVH5_header_t uvh5_header = {0};
      char tel_info_toml_filepath[70] = {'\0'};
      char obs_info_toml_filepath[70] = {'\0'};
      double obs_antenna_positions[BLADE_ATA_INPUT_NANT*3] = {0}, obs_beam_coordinates[BLADE_ATA_OUTPUT_NBEAM*2] = {0};
      double obs_phase_center[2] = {0};
      struct blade_ata_observation_meta observationMetaData = {0};
      struct LonLatAlt arrayReferencePosition = {0};

      double _Complex* antenna_calibration_coeffs;
      char obs_antenna_calibration_filepath[70] = {'\0'};
      char** obs_antenna_names = NULL;

      int fenchan;
      hgetu8(databuf_header, "SCHAN", &observationMetaData.frequencyStartIndex);
      hgetr8(databuf_header, "CHAN_BW", &observationMetaData.channelBandwidthHz);
      hgeti4(databuf_header, "FENCHAN", &fenchan);
      hgetr8(databuf_header, "OBSFREQ", &observationMetaData.rfFrequencyHz);

      double tmp = (double)observationMetaData.rfFrequencyHz +
        (-(double)observationMetaData.frequencyStartIndex - ((double)BLADE_ATA_CONFIG.inputDims.NCHANS / 2.0)
          + ((double)fenchan / 2.0) + 0.5
        ) * (double)observationMetaData.channelBandwidthHz;


      observationMetaData.rfFrequencyHz = tmp;

      observationMetaData.rfFrequencyHz *= 1e6;
      observationMetaData.channelBandwidthHz *= 1e6;
      observationMetaData.totalBandwidthHz = fenchan * observationMetaData.channelBandwidthHz;

      hashpipe_status_lock_safe(user_data->status);
      {
        hgets(user_data->status->buf, "TELINFOP", 70, tel_info_toml_filepath);
        hgets(user_data->status->buf, "OBSINFOP", 70, obs_info_toml_filepath);
        hgets(user_data->status->buf, "CALWGHTP", 70, obs_antenna_calibration_filepath);
      }
      hashpipe_status_unlock_safe(user_data->status);

      hashpipe_info(user_data->thread_name, "Parsing '%s' as Telescope information.", tel_info_toml_filepath);
      UVH5toml_parse_telescope_info(tel_info_toml_filepath, &uvh5_header);
      hashpipe_info(user_data->thread_name, "Parsing '%s' as Observation information.", obs_info_toml_filepath);
      UVH5toml_parse_observation_info(obs_info_toml_filepath, &uvh5_header);
      UVH5Hadmin(&uvh5_header);
      arrayReferencePosition.LAT = calc_rad_from_degree(uvh5_header.latitude);
      arrayReferencePosition.LON = calc_rad_from_degree(uvh5_header.longitude);
      arrayReferencePosition.ALT = uvh5_header.altitude;

      obs_antenna_names = malloc(uvh5_header.Nants_data*sizeof(char*));

      // At this point we have XYZ uvh5_header.antenna_positions, and ENU uvh5_header._antenna_enu_positions
      for(int i = 0; i < uvh5_header.Nants_data; i++) {
        int ant_idx = uvh5_header._antenna_num_idx_map[
          uvh5_header._antenna_numbers_data[i]
        ];
        obs_antenna_positions[i*3 + 0] = uvh5_header.antenna_positions[ant_idx*3 + 0];
        obs_antenna_positions[i*3 + 1] = uvh5_header.antenna_positions[ant_idx*3 + 1];
        obs_antenna_positions[i*3 + 2] = uvh5_header.antenna_positions[ant_idx*3 + 2];

        obs_antenna_names[i] = uvh5_header.antenna_names[ant_idx];
      }
      // BLADE needs ECEF coordinates
      // It doesn't make sense to convert from ECEF back to XYZ in the uvh5 library
      // but then back to ECEF here. TODO don't do the above
      calc_position_to_ecef_frame_from_xyz(
        obs_antenna_positions,
        uvh5_header.Nants_data,
        arrayReferencePosition.LON,
        arrayReferencePosition.LAT,
        arrayReferencePosition.ALT
      );
      // observationMetaData.referenceAntennaIndex = uvh5_header._antenna_num_idx_map[
      //   uvh5_header._antenna_numbers_data[0]
      // ];
      observationMetaData.referenceAntennaIndex = 0;

      collect_beamCoordinates(BLADE_ATA_OUTPUT_NBEAM,
          obs_beam_coordinates, obs_phase_center, databuf_header);
      hashpipe_info(user_data->thread_name, "Parsing '%s' for antenna-weights information.", obs_antenna_calibration_filepath);
      if (
        read_antenna_weights(
          obs_antenna_calibration_filepath,
          uvh5_header.Nants_data, // number of antenna of interest
          obs_antenna_names, // the antenna of interest
          observationMetaData.frequencyStartIndex, // the first channel
          BLADE_ATA_CONFIG.inputDims.NCHANS, // the number of channels
          &antenna_calibration_coeffs // return value
        )
      ) {
        // Failed to open CALWGHTP file, set 1.0+0.0j
        hashpipe_warn(user_data->thread_name, "CALWGHTP `%s` could not be opened. Using 1.0+0.0j.", obs_antenna_calibration_filepath);
        errno = 0;
        size_t size_of_calib =
          BLADE_ATA_CONFIG.inputDims.NANTS*
          BLADE_ATA_CONFIG.inputDims.NCHANS*
          BLADE_ATA_CONFIG.inputDims.NPOLS;
        antenna_calibration_coeffs = malloc(size_of_calib*sizeof(double _Complex*));

        for(int i = 0; i < size_of_calib; i++) {
          antenna_calibration_coeffs[i] = 1.0 + 0.0*I;
        }
      }

      blade_ata_terminate();
      blade_ata_initialize(
        BLADE_ATA_CONFIG,
        user_data->blade_number_of_workers,
        &observationMetaData,
        &arrayReferencePosition,
        obs_phase_center,
        obs_beam_coordinates,
        obs_antenna_positions,
        antenna_calibration_coeffs
      );

      free(obs_antenna_names);
      free(antenna_calibration_coeffs);
    }

    user_data->prev_pktidx_obs_start = pktidx_obs_start;

    if (BLADE_BLOCK_DATA_SIZE != blade_ata_get_output_size()*BLADE_ATA_OUTPUT_ELEMENT_BYTES) {
      hashpipe_error(user_data->thread_name, "BLADE_BLOCK_DATA_SIZE %lu != %lu BLADE configured output size.", BLADE_BLOCK_DATA_SIZE, blade_ata_get_output_size()*BLADE_ATA_OUTPUT_ELEMENT_BYTES);
      pthread_exit(NULL);
      return false;
    }
  }

  return true;
}

bool blade_cb_input_buffer_fetch(void* user_data_void, void** buffer, size_t* buffer_id) {
  blade_userdata_t* user_data = (blade_userdata_t*) user_data_void;

  // prefetch callback has ascertained that the current buffer is filled... don't check again

  #if BLADE_ATA_MODE == BLADE_ATA_MODE_H
  user_data->accumulator_counter = blade_ata_h_accumulator_counter();
  #endif

  *buffer = hpguppi_databuf_data(user_data->in, user_data->in_index);
  clock_gettime(CLOCK_MONOTONIC, &user_data->ts_blocks_recvd[user_data->in_index]);

  blade_ata_set_block_time_mjd(jd_mid_block(hpguppi_databuf_header(user_data->in, user_data->in_index)));
  double dut1;
  hgetr8(hpguppi_databuf_header(user_data->in, user_data->in_index), "DUT1", &dut1);
  blade_ata_set_block_dut1(dut1);

  // hashpipe_info(user_data->thread_name, "batched input block #%d.", user_data->in_index);
  *buffer_id = user_data->in_index;
  user_data->in_index  = (user_data->in_index + 1) % user_data->in->header.n_block;

  return true;
}

void blade_cb_input_buffer_enqueued(void* user_data_void, size_t buffer_id) {
  blade_userdata_t* user_data = (blade_userdata_t*) user_data_void;

  double tbin;
  double obsbw;

  {// Asynchronous CPU work
    // copy across the header
    char* databuf_header = hpguppi_databuf_header(user_data->out, user_data->out_index_free);

    #if BLADE_ATA_MODE == BLADE_ATA_MODE_H
      size_t block_stop_pktidx;
      if (user_data->accumulator_counter == 0) {
        memcpy(databuf_header,
              hpguppi_databuf_header(user_data->in, buffer_id),
              BLOCK_HDR_SIZE);
      }
      else {
        hgetu8(hpguppi_databuf_header(user_data->in, buffer_id), "BLKSTOP", &block_stop_pktidx);
        hputu8(databuf_header, "BLKSTOP", block_stop_pktidx);
      }
    #else
      memcpy(databuf_header,
            hpguppi_databuf_header(user_data->in, buffer_id),
            BLOCK_HDR_SIZE);
    #endif

    //TODO upate output_buffer headers to reflect that they contain beams
    hputi4(databuf_header, "INCOBEAM", (BLADE_ATA_OUTPUT_INCOHERENT_BEAM ? 1 : 0));
    hputi4(databuf_header, "NBEAM", BLADE_ATA_CONFIG.beamformerBeams + (BLADE_ATA_OUTPUT_INCOHERENT_BEAM ? 1 : 0));
    hputi4(databuf_header, "NBITS", BLADE_ATA_OUTPUT_NBITS);
    hputs(databuf_header, "DATATYPE", "FLOAT");
    hputs(databuf_header, "SMPLTYPE", BLADE_ATA_OUTPUT_SAMPLE_TYPE);
    hputi4(databuf_header, "BLOCSIZE", BLADE_BLOCK_DATA_SIZE);

    hgetr8(databuf_header, "TBIN", &tbin);
    hgetr8(databuf_header, "OBSBW", &obsbw);

    #if BLADE_ATA_MODE == BLADE_ATA_MODE_A
    // offload to the downstream filbank writer, which splits OBSNCHAN by number of beams...
    hputi4(databuf_header, "OBSNCHAN", BLADE_ATA_CONFIG.inputDims.NCHANS*BLADE_ATA_CONFIG.channelizerRate*(BLADE_ATA_CONFIG.beamformerBeams + (BLADE_ATA_OUTPUT_INCOHERENT_BEAM ? 1 : 0)));
    hputi4(databuf_header, "NPOL", BLADE_ATA_CONFIG.numberOfOutputPolarizations);

    tbin *= BLADE_ATA_CONFIG.channelizerRate;
    tbin *= BLADE_ATA_CONFIG.integrationSize;

    // negate OBSBW to indicate descending frequency-channel order
    obsbw *= -1.0;

    #elif BLADE_ATA_MODE == BLADE_ATA_MODE_H
    if (user_data->accumulator_counter == 0) {
      // offload to the downstream filbank writer, which splits OBSNCHAN by number of beams...
      hputi4(databuf_header, "OBSNCHAN", BLADE_ATA_CONFIG.inputDims.NCHANS*BLADE_ATA_CONFIG.channelizerRate*BLADE_ATA_CONFIG.accumulateRate*BLADE_ATA_CONFIG.inputDims.NTIME*(BLADE_ATA_CONFIG.beamformerBeams + (BLADE_ATA_OUTPUT_INCOHERENT_BEAM ? 1 : 0)));
      hputi4(databuf_header, "NTIME", 1); // accumulation limits output to NTIME dimension-length of 1
      hputi4(databuf_header, "NPOL", BLADE_ATA_CONFIG.numberOfOutputPolarizations);

      tbin *= BLADE_ATA_CONFIG.channelizerRate;
      tbin *= BLADE_ATA_CONFIG.accumulateRate*BLADE_ATA_CONFIG.inputDims.NTIME;

      // negate OBSBW to indicate descending frequency-channel order
      obsbw *= -1.0;
    }

    #else
    hputi4(databuf_header, "NCHAN", BLADE_ATA_CONFIG.inputDims.NCHANS*BLADE_ATA_CONFIG.channelizerRate); // beams are split into separate files...
    hputi4(databuf_header, "OBSNCHAN", BLADE_ATA_CONFIG.inputDims.NCHANS*BLADE_ATA_CONFIG.channelizerRate); // beams are split into separate files...
    #endif

    hputr8(databuf_header, "TBIN", tbin);
    hputr8(databuf_header, "OBSBW", obsbw);
  }

}

void blade_cb_input_buffer_ready(void* user_data_void, const void* buffer, size_t buffer_id) {
  blade_userdata_t* user_data = (blade_userdata_t*) user_data_void;

  hpguppi_databuf_set_free(user_data->in, buffer_id);
  // hashpipe_info(user_data->thread_name, "freed inputput block #%d.", buffer_id);
}

bool blade_cb_output_buffer_fetch(void* user_data_void, void** buffer, size_t* buffer_id) {
  blade_userdata_t* user_data = (blade_userdata_t*) user_data_void;

  {// poll if output buffer is free
    int hpguppi_databuf_wait_rv = hpguppi_databuf_wait_free_timeout(
      user_data->out, user_data->out_index_free,
      &user_data->ts_buffer_wait_timeout
    );
    if (hpguppi_databuf_wait_rv == HASHPIPE_TIMEOUT) {
      if (user_data->status_state != 2)
      {
        hashpipe_status_lock_safe(user_data->status);
        {
          hputs(user_data->status->buf, user_data->status_key, "outblocked");
        }
        hashpipe_status_unlock_safe(user_data->status);
        user_data->status_state = 2;
      }
      return false;
    }
    if(hpguppi_databuf_wait_rv != HASHPIPE_OK) {
      hashpipe_error(user_data->thread_name, "error waiting for output buffer #%d, rv: %i", user_data->out_index_free, hpguppi_databuf_wait_rv);
      pthread_exit(NULL);
      return false;
    }
  }

  #if BLADE_ATA_MODE == BLADE_ATA_MODE_A || BLADE_ATA_MODE == BLADE_ATA_MODE_H
  *buffer = user_data->out_intermediary[user_data->out_index_free];
  #else
  hpguppi_blade_output_databuf_t* out = (hpguppi_blade_output_databuf_t*) user_data->out;
  *buffer = hpguppi_databuf_data(out, user_data->out_index_free);
  #endif

  // hashpipe_info(user_data->thread_name, "batched output block #%d.", user_data->out_index_free);
  *buffer_id = user_data->out_index_free;
  user_data->out_index_free  = (user_data->out_index_free + 1) % user_data->out->header.n_block;
  return true;
}

void blade_cb_output_buffer_ready(void* user_data_void, const void* buffer, size_t buffer_id) {
  blade_userdata_t* user_data = (blade_userdata_t*) user_data_void;

  #if BLADE_ATA_MODE == BLADE_ATA_MODE_A || BLADE_ATA_MODE == BLADE_ATA_MODE_H
  // Mode A/H has filterbank outputs, so transpose the binary data to suit that of .fil files
  const int npol = BLADE_ATA_CONFIG.numberOfOutputPolarizations;
  const int nbeams = BLADE_ATA_CONFIG.beamformerBeams + (BLADE_ATA_OUTPUT_INCOHERENT_BEAM ? 1 : 0);

  #if BLADE_ATA_MODE == BLADE_ATA_MODE_A
  float* ftp_output = user_data->out_intermediary[buffer_id];
  float* tpf_output = (float*) hpguppi_databuf_data(user_data->out, buffer_id);
  const int nfreq = BLADE_ATA_CONFIG.inputDims.NCHANS*BLADE_ATA_CONFIG.channelizerRate;
  const int ntime = BLADE_ATA_CONFIG.inputDims.NTIME / (BLADE_ATA_CONFIG.integrationSize * BLADE_ATA_CONFIG.channelizerRate);
  #elif BLADE_ATA_MODE == BLADE_ATA_MODE_H
  float* ftp_output = user_data->out_intermediary[buffer_id];
  float* tpf_output = (float*) hpguppi_databuf_data(user_data->out, buffer_id);
  const int nfreq = BLADE_ATA_CONFIG.inputDims.NCHANS*BLADE_ATA_CONFIG.channelizerRate*BLADE_ATA_CONFIG.accumulateRate*BLADE_ATA_CONFIG.inputDims.NTIME;
  const int ntime = 1;
  #endif
  int b,f,t,p;

  for(b = 0; b < nbeams; b++) {
    for(f = 0; f < nfreq; f++) {
      for(t = 0; t < ntime; t++) {
        for(p = 0; p < npol; p++) {
          tpf_output[
              ((  b *ntime
                + t)*npol
                + p)*nfreq
                - f + nfreq-1 // filterbank files typically are frequency descending
          ] = ftp_output[
              ((  b *nfreq
                + f)*ntime
                + t)*npol
                + p
          ];
        }
      }
    }
  }
  #endif

  // hashpipe_info(user_data->thread_name, "filled output block #%d.", buffer_id);
  if(user_data->out_index_fill != buffer_id) {
    hashpipe_error(user_data->thread_name, "Non sequential output buffer fill detected: %d != %d.", user_data->out_index_fill, buffer_id);
    pthread_exit(NULL);
    return;
  }
  hpguppi_databuf_set_filled(user_data->out, buffer_id);
  user_data->out_index_fill = (user_data->out_index_fill + 1) % user_data->out->header.n_block;

  // Update moving sum (for moving average)
  struct timespec ts_free_input = {0};
  clock_gettime(CLOCK_MONOTONIC, &ts_free_input);

  uint64_t fill_to_free_elapsed_ns = ELAPSED_NS(user_data->ts_blocks_recvd[buffer_id], ts_free_input);
  user_data->fill_to_free_moving_sum_ns +=
      fill_to_free_elapsed_ns - user_data->fill_to_free_block_ns[buffer_id];
  // Store new value
  user_data->fill_to_free_block_ns[buffer_id] = fill_to_free_elapsed_ns;
}

void blade_cb_clear_queued_input(void* user_data_void, size_t input_id) {
  blade_userdata_t* user_data = (blade_userdata_t*) user_data_void;
  // hashpipe_info(user_data->thread_name, "cancelled input block #%d.", input_id);
  hpguppi_databuf_set_free(user_data->in, input_id);
}

void blade_cb_clear_queued_output(void* user_data_void, size_t output_id) {
  blade_userdata_t* user_data = (blade_userdata_t*) user_data_void;
  if(user_data->out_index_free != user_data->out_index_fill) {
    hashpipe_info(user_data->thread_name, "reset output block index from #%d to #%d.", user_data->out_index_free, user_data->out_index_fill);
    user_data->out_index_free = user_data->out_index_fill;
  }
}

static void *run(hashpipe_thread_args_t *args)
{
  blade_userdata_t blade_userdata = {
    .blade_number_of_workers = 2,

    .status = &args->st,

    .thread_name = args->thread_desc->name,
    .status_key = args->thread_desc->skey,

    .status_state = 0,
    .ts_last_status_update = {
      .tv_sec = 0,
      .tv_nsec = 0,
    },
    .ts_buffer_wait_timeout = {
      .tv_sec = 0,
      .tv_nsec = 500000000, // 500 ms
    },


    .prev_flagged_NANTS = 0,
    .prev_flagged_NCHAN = 0,
    .prev_flagged_NTIME = 0,
    .prev_flagged_NPOLS = 0,

    .accumulator_counter = 0,

    .in_index = 0,
    .out_index_free = 0,

    .in = (hpguppi_input_databuf_t *)args->ibuf,
    .out = (hpguppi_blade_output_databuf_t *)args->obuf,

    // .out_intermediary = {NULL},

    .fill_to_free_moving_sum_ns = 0,
    // .fill_to_free_block_ns = {0},
    // .ts_blocks_recvd = {0},
  };

  #if BLADE_ATA_MODE == BLADE_ATA_MODE_A || BLADE_ATA_MODE == BLADE_ATA_MODE_H
  for(size_t i = 0; i < N_INPUT_BLOCKS; i++) {
    blade_userdata.out_intermediary[i] = malloc(BLADE_BLOCK_OUTPUT_DATA_SIZE);
    blade_userdata.fill_to_free_block_ns[i] = 0;
    memset(blade_userdata.ts_blocks_recvd + i, 0, sizeof(struct timespec));
  }
  #endif

  int cudaDeviceId = args->instance_id;

  hashpipe_status_lock_safe(blade_userdata.status);
  {
    hgeti4(blade_userdata.status->buf, "CUDADEV", &cudaDeviceId);
  }
  hashpipe_status_unlock_safe(blade_userdata.status);

  if (cudaDeviceId >= 0) {
    if (blade_use_device(cudaDeviceId)) {
      hashpipe_info(args->thread_desc->name, "Successfully set CUDA device to %d.", cudaDeviceId);
    }
    else {
      hashpipe_info(args->thread_desc->name, "Failed to set CUDA device to %d.", cudaDeviceId);
      cudaDeviceId = -1;
    }
  }
  hashpipe_status_lock_safe(blade_userdata.status);
  {
    hputi4(blade_userdata.status->buf, "CUDADEV", cudaDeviceId);
    hputi4(blade_userdata.status->buf, "BLDBLKSZ", BLADE_BLOCK_DATA_SIZE);
  }
  hashpipe_status_unlock_safe(blade_userdata.status);

  for(int i = 0; i < N_INPUT_BLOCKS; i++)
  {
    blade_pin_memory(hpguppi_databuf_data(blade_userdata.in, i), BLOCK_DATA_SIZE);
    blade_pin_memory(hpguppi_databuf_data(blade_userdata.out, i), BLADE_BLOCK_DATA_SIZE);
  }

  blade_ata_register_user_data((void*) &blade_userdata);
  blade_ata_register_input_buffer_prefetch_cb(&blade_cb_input_buffer_prefetch);
  blade_ata_register_input_buffer_fetch_cb(&blade_cb_input_buffer_fetch);
  blade_ata_register_input_buffer_ready_cb(&blade_cb_input_buffer_ready);
  blade_ata_register_input_buffer_enqueued_cb(&blade_cb_input_buffer_enqueued);
  blade_ata_register_output_buffer_fetch_cb(&blade_cb_output_buffer_fetch);
  blade_ata_register_output_buffer_ready_cb(&blade_cb_output_buffer_ready);
  blade_ata_register_blade_queued_input_clear_cb(&blade_cb_clear_queued_input);
  blade_ata_register_blade_queued_output_clear_cb(&blade_cb_clear_queued_output);

  while (run_threads())
  {
    blade_ata_compute_step();

    // Will exit if thread has been cancelled
    pthread_testcancel();
  }

  hashpipe_info(blade_userdata.thread_name, "returning");
  blade_ata_terminate();
  return NULL;
}

static hashpipe_thread_desc_t blade_beamformer_thread = {
  name: "hpguppi_ata_blade_beamformer_thread",
  skey: "BEAMSTAT",
  init: NULL,
  run: run,
  ibuf_desc: {hpguppi_input_databuf_create},
  obuf_desc: {hpguppi_blade_output_databuf_create}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(& blade_beamformer_thread);
}
