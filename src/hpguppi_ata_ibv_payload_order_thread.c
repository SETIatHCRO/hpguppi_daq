// hpguppi_atasnap_ibv_payload_order_thread.c
//
// A Hashpipe thread that proceeses "voltage mode" packets sent by the ATA SNAP
// design from an input buffer (populated by hpguppi_ibverbs_pkt_thread) and
// assembles them into GUPPI RAW blocks.

// TODO TEST Wait for first (second?) start-of-block when transitioning into
//           LISTEN state so that the first block will be complete.
// TODO Add PSPKTS and PSDRPS status buffer fields for pktsock
// TODO TEST Set NETSTAE to idle in IDLE state
// TODO TEST IP_DROP_MEMBERSHIP needs mcast IP address (i.e. not 0.0.0.0)

#define _GNU_SOURCE 1
//#include <stdio.h>
//#include <sys/types.h>
#include <stdlib.h>
#include <sched.h>
#include <math.h>
#include <unistd.h>
#include <limits.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <poll.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "hashpipe.h"
#include "hpguppi_databuf.h"
#include "hpguppi_time.h"
#include "hpguppi_util.h"
#include "hpguppi_atasnap.h"
#include "hpguppi_pktbuf.h"

#include <omp.h>
#define ATA_IBV_FOR_PACKET_THREAD_COUNT 10
#define ATA_IBV_TRANSPOSE_PACKET_THREAD_COUNT 1
#define ATA_IBV_THREAD_COUNT ATA_IBV_FOR_PACKET_THREAD_COUNT*ATA_IBV_TRANSPOSE_PACKET_THREAD_COUNT

// Change to 1 to use temporal memset() rather than non-temporal bzero_nt()
#if 0
#define bzero_nt(d,l) memset(d,0,l)
#endif

// Change to 1 to use temporal memcpy() rather than non-temporal memcpy_nt()
#if 0
#define memcpy_nt(dst,src,len) memcpy(dst,src,len)
#endif

// This thread's init() function, if provided, is called by the Hashpipe
// framework at startup to allow the thread to perform initialization tasks
// such as setting up network connections or GPU devices.
static int init(hashpipe_thread_args_t *args)
{
  // Local aliases to shorten access to args fields
  // Our input buffer happens to be a hpguppi_input_databuf
  hpguppi_input_databuf_t *dbin  = (hpguppi_input_databuf_t *)args->ibuf;
  const char * thread_name = args->thread_desc->name;
  const char * status_key = args->thread_desc->skey;
  hashpipe_status_t *st = &args->st;

  // THREAD_COUNT align number of packet-slots per block (ibverbs_pkt_thread utilised)
  struct hpguppi_pktbuf_info * pktbuf_info = hpguppi_pktbuf_info_ptr(dbin);
  hashpipe_info(thread_name, "Aligning %d packet-slots per block to nearest smaller multiple of %d", pktbuf_info->slots_per_block, ATA_IBV_THREAD_COUNT);
  pktbuf_info->slots_per_block -= pktbuf_info->slots_per_block%ATA_IBV_THREAD_COUNT;

  // Non-network essential paramaters
  int blocsize=BLOCK_DATA_SIZE;
  int directio=1;
  int nbits=4;
  int npol=4;
  double obsfreq=0;
  double chan_bw=900.0/4096;
  double obsbw=256*chan_bw;
  int obsnchan=1;
  int nants=1;
  int overlap=0;
  double tbin=1e-6;
  char obs_mode[80] = {0};
  struct rlimit rlim;
  // Port to listen on
  uint32_t bindport = 0; // zero leads to no sniffer_flow being created

  strcpy(obs_mode, "RAW");

  // Verify that the IBVPKTSZ was specified as expected/requried
  if(hpguppi_pktbuf_slot_offset(dbin, ATA_SNAP_PKT_OFFSET_HEADER) %
      PKT_ALIGNMENT_SIZE != 0
  || hpguppi_pktbuf_slot_offset(dbin, ATA_SNAP_PKT_OFFSET_PAYLOAD) %
      PKT_ALIGNMENT_SIZE != 0) {
    errno = EINVAL;
    hashpipe_error(thread_name, "IBVPKTSZ!=%d,%d,[...]",
        ATA_SNAP_PKT_OFFSET_HEADER, ATA_SNAP_PKT_SIZE_HEADER);
    return HASHPIPE_ERR_PARAM;
  }

  // Set RLIMIT_RTPRIO to 1
  getrlimit(RLIMIT_RTPRIO, &rlim);
	rlim.rlim_cur = 1;
	rlim.rlim_max = (rlim.rlim_max < rlim.rlim_cur ? rlim.rlim_cur : rlim.rlim_max);
  if(setrlimit(RLIMIT_RTPRIO, &rlim)) {
    hashpipe_error(thread_name, "setrlimit(RLIMIT_RTPRIO)");
  }

  struct sched_param sched_param = {
    .sched_priority = 1
  };
  if(sched_setscheduler(0, SCHED_RR, &sched_param)) {
    hashpipe_error(thread_name, "sched_setscheduler");
  }

  hashpipe_status_lock_safe(st);
  {
    // Get info from status buffer if present (no change if not present)
    hgetu4(st->buf, "BINDPORT", &bindport);

    hgeti4(st->buf, "BLOCSIZE", &blocsize);
    hgeti4(st->buf, "DIRECTIO", &directio);
    hgeti4(st->buf, "NANTS", &nants);
    hgeti4(st->buf, "NBITS", &nbits);
    hgeti4(st->buf, "NPOL", &npol);
    hgetr8(st->buf, "OBSFREQ", &obsfreq);
    hgetr8(st->buf, "OBSBW", &obsbw);
    hgetr8(st->buf, "CHAN_BW", &chan_bw);
    hgeti4(st->buf, "OBSNCHAN", &obsnchan);
    hgeti4(st->buf, "OVERLAP", &overlap);
    hgets(st->buf, "OBS_MODE", sizeof(obs_mode), obs_mode);

    // Prevent div-by-zero errors (should never happen...)
    if(nants == 0) {
      nants = 1;
    }

    // If CHAN_BW is zero, set to default value (1 MHz)
    if(chan_bw == 0.0) {
      chan_bw = 1.0;
    }

    // Calculate tbin and obsbw from chan_bw
    tbin = 1e-6 / fabs(chan_bw);
    obsbw = chan_bw * obsnchan / nants;

    // Update status buffer (in case fields were not there before).
    hputs(st->buf, "DAQSTATE", "LISTEN");
    hputi4(st->buf, "BLOCSIZE", blocsize);
    hputi4(st->buf, "DIRECTIO", directio);
    hputi4(st->buf, "NANTS", nants);
    hputi4(st->buf, "NBITS", nbits);
    hputi4(st->buf, "NPOL", npol);
    hputr8(st->buf, "OBSBW", obsbw);
    hputr8(st->buf, "CHAN_BW", chan_bw);
    hputi4(st->buf, "OBSNCHAN", obsnchan);
    hputi4(st->buf, "OVERLAP", overlap);
    hputs(st->buf, "PKTFMT", "ATASNAPV");
    hputr8(st->buf, "TBIN", tbin);
    hputs(st->buf, "OBS_MODE", obs_mode);
    hputi4(st->buf, "NDROP", 0);
    // Set status_key to init
    hputs(st->buf, status_key, "init");
    hputu8(st->buf, "PKTSTART", 0);
    hputu8(st->buf, "PKTSTOP", 0);

    // Store bind port in status buffer (in case it was not there before).
    hputu4(st->buf, "BINDPORT", bindport);
    
    if(bindport > 0) {
      hashpipe_info(thread_name, "Setting IBVSNIFF to BINDPORT (%d) to trigger sniffer ibv_flow creation.", bindport);
      hputi4(st->buf, "IBVSNIFF", bindport);
    }
  }
  hashpipe_status_unlock_safe(st);

  // Success!
  return 0;
}

static void * run(hashpipe_thread_args_t * args)
{
#if 0
int debug_i=0, debug_j=0;
#endif
  // Local aliases to shorten access to args fields
  // Our input and output buffers happen to be a hpguppi_input_databuf
  hpguppi_input_databuf_t *dbin  = (hpguppi_input_databuf_t *)args->ibuf;
  ATA_IBV_OUT_DATABUF_T *dbout = (ATA_IBV_OUT_DATABUF_T *)args->obuf;
  hashpipe_status_t *st = &args->st;
  const char * thread_name = args->thread_desc->name;
  const char * status_key = args->thread_desc->skey;

  // Update status_key with idle state
  hashpipe_status_lock_safe(st);
  {
    hputs(st->buf, status_key, "listen");
  }
  hashpipe_status_unlock_safe(st);

  // Misc counters, etc
  int rv=0;
  int i;//, pkt_chan_idx;

#if 0
  uint64_t u64;
  uint8_t u8 = 0;
  uint8_t *pu8in = (uint8_t *)dbin;
  uint8_t *pu8out = (uint8_t *)dbout;
  for(u64=0; u64<sizeof(hpguppi_input_databuf_t); u64+=4096) {
    if(u8 || !u8) {
      u8 += pu8in[u64];
      u8 += pu8out[u64];
    }
  }
  hashpipe_info(thread_name, "db pagein sum is %u", u8);
#endif
  memset(dbout->block, 0, sizeof(dbout->block));
  hashpipe_info(thread_name,
      "set %lu bytes in dbout to 0", sizeof(dbout->block));

  // for(i=0; i<N_INPUT_BLOCKS; i++) {
  //   hashpipe_info(thread_name, "db_in  block %2d : %p %p", i,
  //       hpguppi_databuf_data(dbin, i),
  //       hpguppi_databuf_data(dbin, i) + BLOCK_DATA_SIZE - 1);
  // }

  // for(i=0; i<N_INPUT_BLOCKS; i++) {
  //   hashpipe_info(thread_name, "db_out block %2d : %p %p", i,
  //       hpguppi_databuf_data(dbout, i),
  //       hpguppi_databuf_data(dbout, i) + BLOCK_DATA_SIZE - 1);
  // }

  // The incoming packets are taken from blocks of the input databuf and then
  // converted to GUPPI RAW format in blocks of the output databuf to pass to
  // the downstream thread.  We currently support two active output blocks (aka
  // "working blocks").  Working blocks are associated with absolute output
  // block numbers, which are simply PKTIDX values divided by the number of
  // packets per block (discarding any remainder).  Let the block numbers for
  // the first working block (wblk[0]) be W.  The block number for the second
  // working block (wblk[1]) will be W+1.  Incoming packets corresponding to
  // block W or W+1 are placed in the corresponding data buffer block.
  // Incoming packets for block W+2 cause block W to be "finalized" and handed
  // off to the downstream thread, working block 1 moves to working block 0 and
  // working block 1 is incremented to be W+2.  Things get "interesting" when a
  // packet is recevied for block < W or block > W+2.  Packets for block W-1
  // are ignored.  Packets with PKTIDX P corresponding block < W-1 or block >
  // W+2 cause the current working blocks' block numbers to be reset such that
  // W will refer to the block containing P and W+1 will refer to the block
  // after that.
  //
  // wblk is a two element array of block_info structures (i.e. the working
  // blocks)
  int wblk_idx;
  const int n_wblock = 3;
  struct datablock_stats wblk[n_wblock];
  uint32_t *thread_wblk_pkt_count = malloc(ATA_IBV_FOR_PACKET_THREAD_COUNT*n_wblock*sizeof(uint32_t));
  memset(thread_wblk_pkt_count, 0, ATA_IBV_FOR_PACKET_THREAD_COUNT*n_wblock*sizeof(uint32_t));

  // Packet block variables
  uint64_t blk0_relative_pkt_seq_num = 0;
  unsigned long pkt_blk_num, last_pkt_blk_num = ~0;
  uint64_t obs_start_seq_num = 0, obs_stop_seq_num = 0, blk0_start_seq_num = 0;
  uint64_t prev_obs_start_seq_num, prev_obs_stop_seq_num;
  uint32_t antenna_byte_stride, channel_byte_stride, time_byte_stride;

  // Heartbeat variables
  struct timespec ts_start_block = {0}, ts_stop_block = {0};
  struct timespec ts_checked_obs_info = {0}, ts_tried_obs_info = {0}, ts_now = {0};
  const uint64_t obs_info_refresh_period_ns = 200*1000*1000;
  const uint64_t obs_info_retry_period_s = 5;
  uint64_t obs_info_refresh_elapsed_ns;
  float blocks_per_second = 0.0;

  // State variables
  enum obs_info_validity obs_info_validity = OBS_UNKNOWN;
  char buf_status[80];
  int buf_full;  
  
  char waiting=-1, observing=0;
  char flag_reinit_blks=0;
  int32_t observation_complete=0;

  char flag_state_update = 0;
  char  LATE_PKTIDX_flagged = 0;
  char  PKT_OBS_FENG_flagged,
        PKT_OBS_SCHAN_flagged,
        PKT_OBS_NCHAN_flagged,
        PKT_OBS_PKTNTIME_flagged,
        PKT_OBS_PKTIDX_flagged;

  // Variables for working with the input databuf
  struct hpguppi_pktbuf_info * pktbuf_info = hpguppi_pktbuf_info_ptr(dbin);
  int block_idx_in = 0;
  struct timespec timeout_in = {0, 50 * 1000 * 1000}; // 50 ms
  char * datablock_header;

  // Variables for counting packets and bytes.
  uint64_t npacket=0, npacket_total=0, npacket_drop=0, ndrop_total=0;
  uint64_t n_blks_rushed=0;
  // uint32_t nbogus_size;

  // Variables for handing received packets
  uint8_t * p_u8pkt;
  struct ata_snap_ibv_pkt * p_pkt = NULL;
  PKT_PAYLOAD_CP_T* p_payload = NULL;
  PKT_PAYLOAD_CP_T* dest_feng_pktidx_offset = NULL;

  // Structure to hold observation info, init all fields to invalid values
  struct ata_snap_obs_info obs_info;
  ata_snap_obs_info_init(&obs_info);

  // Structure to hold feng info from packet
  struct ata_snap_pkt_info pkt_info = {0};

  // Variables for tracking timing stats
  //
  // ts_start_recv(N) to ts_stop_recv(N) is the time spent in the "receive" call.
  // ts_stop_recv(N) to ts_start_recv(N+1) is the time spent processing received data.
  struct timespec ts_start_recv = {0}, ts_stop_recv = {0};
  struct timespec ts_free_input = {0};

  // Used to calculate moving average of fill-to-free times for input blocks
  uint64_t fill_to_free_elapsed_ns;
  uint64_t fill_to_free_moving_sum_ns = 0;
  uint64_t fill_to_free_block_ns[N_INPUT_BLOCKS] = {0};

  //struct timespec ts_sleep = {0, 10 * 1000 * 1000}; // 10 ms

#if 0
  // Allocate a 2K buffer into which packet will be non-temporally copied
  // before processing.  This buffer will be cached (due to parsing of the
  // headers), but the input databuf blocks will not be cached.
  if((rv = posix_memalign((void **)&p_spdpkt, 4096, MAX_PKT_SIZE))) {
    errno = rv;
    hashpipe_error(thread_name, "cannot allocate page aligned packet buffer");
    return NULL;
  }
#endif

  // Initialize working blocks
  for(wblk_idx=0; wblk_idx<n_wblock; wblk_idx++) {
      // important to initialise the blocks with unique block idx
      // and to set the block_num < -1 so the first observation packet causes a 
      // discontinuity and so a contemporary re-init of the working blocks
      init_datablock_stats(wblk+wblk_idx, dbout, wblk_idx, -10, obs_info.pkt_per_block);
      wait_for_block_free(wblk+wblk_idx, st, status_key);
  }

  // Get any obs info from status buffer, store values
  ata_snap_obs_info_read(st, &obs_info);
  ata_snap_obs_info_write(st, &obs_info);

  // Wait for ibvpkt thread to be running, then it's OK to add/remove flows.
  hpguppi_ibvpkt_wait_running(st);

  const size_t slots_per_block = pktbuf_info->slots_per_block;
  const size_t slot_size = pktbuf_info->slot_size;
  if(slots_per_block % ATA_IBV_THREAD_COUNT != 0 ){
    hashpipe_warn(thread_name, "The slots per block (%lu) should be a multiple of the parallelism (%d) for optimal performance.", slots_per_block, ATA_IBV_THREAD_COUNT);
  }

  hashpipe_status_lock_safe(st);
  {
    hputi4(st->buf, "NETTHRDS", ATA_IBV_THREAD_COUNT);
  }
  hashpipe_status_unlock_safe(st);

  #if ATA_IBV_TRANSPOSE_PACKET_THREAD_COUNT > 1 && ATA_IBV_FOR_PACKET_THREAD_COUNT > 1
    omp_set_nested(1);
  #endif

  // set OBSDONE in case no downstream thread controlls it.
  hput_obsdone(st, 1);
  hashpipe_info(thread_name, "PKT_PAYLOAD_CP_T size: %d", sizeof(PKT_PAYLOAD_CP_T));

  // Main loop
  while (run_threads()) {

    // Mark ts_stop_recv as unset
    ts_stop_recv.tv_sec = 0;

    // Wait for data
    do {
      clock_gettime(CLOCK_MONOTONIC, &ts_start_recv);
      // If ts_stop_recv has been set
      // if(ts_stop_recv.tv_sec != 0) {
      //   // Accumulate processing time
      //   ns_processed_net += ELAPSED_NS(ts_stop_recv, ts_start_recv);
      // }
      rv = hpguppi_databuf_wait_filled_timeout(
          dbin, block_idx_in, &timeout_in);
      clock_gettime(CLOCK_MONOTONIC, &ts_stop_recv);
      memcpy(&ts_now, &ts_stop_recv, sizeof(struct timespec));

      obs_info_refresh_elapsed_ns = ELAPSED_NS(ts_checked_obs_info, ts_now);


      if(rv == HASHPIPE_TIMEOUT && obs_info_refresh_elapsed_ns <= obs_info_refresh_period_ns) {
        // No, continue receiving
        continue;
      }

      // Got packets or new second

      // We perform some status buffer updates every second
      if(flag_state_update || obs_info_refresh_elapsed_ns > obs_info_refresh_period_ns) {
        memcpy(&ts_checked_obs_info, &ts_now, sizeof(struct timespec));

        // write obs_info to overwrite any changes
        if (obs_info_validity == OBS_VALID && // if obs_info is valid
            observing ){ //and observing
            ata_snap_obs_info_write_with_validity(st, &obs_info, obs_info_validity);
        }
        else {//otherwise read obs_info
          if(ata_snap_obs_info_read_with_validity(st, &obs_info, &obs_info_validity)){
            // this code executes if the obs_info has CHANGED
            // (ie at least once before valid observation)
            ata_snap_obs_info_write_with_validity(st, &obs_info, obs_info_validity);
            
            #if ATA_PAYLOAD_TRANSPOSE == ATA_PAYLOAD_TRANSPOSE_FTP
            time_byte_stride = ATASNAP_DEFAULT_PKT_SAMPLE_BYTE_STRIDE;
            channel_byte_stride = obs_info.pktidx_per_block*time_byte_stride;
            antenna_byte_stride = obs_info.nchan*channel_byte_stride;
            #elif ATA_PAYLOAD_TRANSPOSE == ATA_PAYLOAD_TRANSPOSE_TFP
            antenna_byte_stride = ATASNAP_DEFAULT_PKT_SAMPLE_BYTE_STRIDE;
            channel_byte_stride = XGPU_BLOCK_NANTS*antenna_byte_stride;
            time_byte_stride = obs_info.nchan*channel_byte_stride;
            #elif ATA_PAYLOAD_TRANSPOSE == ATA_PAYLOAD_TRANSPOSE_TFP_DP4A
            antenna_byte_stride = ATASNAP_DEFAULT_PKTNPOL*4*ATASNAP_DEFAULT_SAMPLE_BYTESIZE; // `*4` keeps databuf offset logic uniform
            channel_byte_stride = XGPU_BLOCK_NANTS*antenna_byte_stride;
            time_byte_stride = obs_info.nchan*channel_byte_stride/4; // `/4` keeps databuf offset logic uniform
            #endif
            antenna_byte_stride /= sizeof(PKT_PAYLOAD_CP_T);
            channel_byte_stride /= sizeof(PKT_PAYLOAD_CP_T);
            time_byte_stride /= sizeof(PKT_PAYLOAD_CP_T);
            memcpy(&ts_tried_obs_info, &ts_now, sizeof(struct timespec));
          }
          else if (ELAPSED_S(ts_tried_obs_info, ts_now) > obs_info_retry_period_s){
            memcpy(&ts_tried_obs_info, &ts_now, sizeof(struct timespec));
            obs_info_validity = OBS_SEEMS_VALID;
            
            hashpipe_status_lock_safe(st);
              hputs(st->buf, "OBSINFO", "VALID");
            hashpipe_status_unlock_safe(st);

            PKT_OBS_FENG_flagged = 0;
            PKT_OBS_SCHAN_flagged = 0;
            PKT_OBS_NCHAN_flagged = 0;
            PKT_OBS_PKTNTIME_flagged = 0;
            PKT_OBS_PKTIDX_flagged = 0;
          }
        }

        // kill any observations if OBS_info_INVALID
        if (obs_info_validity < OBS_SEEMS_VALID && obs_stop_seq_num > 0){
          obs_start_seq_num = 0;
          obs_stop_seq_num = 0;
          hashpipe_status_lock_safe(st);
          {
            hputu8(st->buf, "PKTSTART", 0);
            hputu8(st->buf, "PKTSTOP", 0);
          }
          hashpipe_status_unlock_safe(st);
        }
        else{
          prev_obs_start_seq_num = obs_start_seq_num;
          prev_obs_stop_seq_num = obs_stop_seq_num;
          hashpipe_status_lock_safe(st);
          {
            hgetu8(st->buf, "PKTSTART", &obs_start_seq_num);
            hgetu8(st->buf, "PKTSTOP", &obs_stop_seq_num);
          }
          hashpipe_status_unlock_safe(st);
          
          if(obs_start_seq_num != prev_obs_start_seq_num){
            hashpipe_info(thread_name, "obs_start_seq_num changed %lu -> %lu", prev_obs_start_seq_num, obs_start_seq_num);
            if(observing){
              hashpipe_info(thread_name, "obs_start_seq_num change ignored while in observation.");
              obs_start_seq_num = prev_obs_start_seq_num;
            }
          }
          if(obs_stop_seq_num != prev_obs_stop_seq_num){
            hashpipe_info(thread_name, "obs_stop_seq_num changed %lu -> %lu", prev_obs_stop_seq_num, obs_stop_seq_num);
          }
        }

        npacket_total += npacket;

        hashpipe_status_lock_safe(st);
        {
          hputi8(st->buf, "NPKTS", npacket_total);
          hgetu8(st->buf, "NDROP", &ndrop_total);
          ndrop_total += npacket_drop;
          hputi8(st->buf, "NDROP", ndrop_total);
          hputr4(st->buf, "PHYSPKPS", npacket*(1e9/obs_info_refresh_elapsed_ns));
          hputr4(st->buf, "PHYSGBPS", (npacket*obs_info.pkt_data_size)/((float) obs_info_refresh_elapsed_ns));

          hputu4(st->buf, "RUSHBLKS", n_blks_rushed);

          hputr4(st->buf, "NETBLKPS", blocks_per_second);
          hputr4(st->buf, "NETBLKMS",
              round((double)fill_to_free_moving_sum_ns / N_INPUT_BLOCKS) / 1e6);

          buf_full = hpguppi_databuf_total_status(dbout);
          sprintf(buf_status, "%d/%d", buf_full, dbout->header.n_block);
          hputs(st->buf, "NETBUFST", buf_status);
        }
        hashpipe_status_unlock_safe(st);
        npacket = 0;
        npacket_drop = 0;
      } // curtime != lasttime

      // Set status field to "waiting" if we are not getting packets
      if(rv && run_threads() && !waiting) {
        hashpipe_status_lock_safe(st);
        {
          hputs(st->buf, status_key, "waiting");
        }
        hashpipe_status_unlock_safe(st);
        waiting=1;
      }

      // Will exit if thread has been cancelled
      pthread_testcancel();
    } while (rv && run_threads()); // end wait for data loop

    if(!run_threads()) {
      // We're outta here!
      // But first mark the block free if we got one.
      if(!rv) {
        hpguppi_databuf_set_free(dbin, block_idx_in);
        clock_gettime(CLOCK_MONOTONIC, &ts_free_input);
        fprintf(stderr, "final fill-to-free %ld ns\n", ELAPSED_NS(ts_stop_recv, ts_free_input));
      }
      break;
    }

    // Got packet(s)!  Update status if needed.
    if (waiting) {
      hashpipe_status_lock_safe(st);
      {
        hputs(st->buf, status_key, "receiving");
      }
      hashpipe_status_unlock_safe(st);
      waiting=0;
    }

    // collect downstream controlled flag
    hget_obsdone(st, &observation_complete);

    // Gather pkt_0 pointer
    p_u8pkt = (uint8_t *)hpguppi_databuf_data(dbin, block_idx_in);

    // parse the first packet to determine whether or not to idle later.
    p_pkt = (struct ata_snap_ibv_pkt *)p_u8pkt;
    ata_snap_parse_ibv_packet(p_pkt, &pkt_info);

    // Only do the work if packets seem to be in range, or the downstream controlled `observation_complete` is low
    if((pkt_info.pktidx <= obs_stop_seq_num || !observation_complete) && obs_info_validity > OBS_INVALID) {
      // Check the first packet's timestamp, to determine reinit_blocks
      //  This works because it is figured that a block filled with packets
      //  will contain packets with indices that place them in two adjacent downstream blocks.
      
      // Get packet index and absolute block number for packet
      blk0_relative_pkt_seq_num = pkt_info.pktidx - blk0_start_seq_num;
      // Get packet's block number relative to the first block's starting index.
      pkt_blk_num = blk0_relative_pkt_seq_num / obs_info.pktidx_per_block;

      //TODO dont use pkt_blk_num due to underflow
      if(pkt_blk_num + 1 < wblk[0].block_num 
        || pkt_blk_num > wblk[n_wblock-1].block_num + 1) {
          if(!observing && observation_complete) {
            flag_reinit_blks = 1;
            blk0_start_seq_num = pkt_info.pktidx;
            align_blk0_with_obsstart(&blk0_start_seq_num, obs_start_seq_num, obs_info.pktidx_per_block);

            // Should only happen when seeing first packet when obs_info is valid
            // warn in case it happens in other scenarios
            hashpipe_warn(thread_name,
                "working blocks reinit due to packet index out of working range\n\t\t(PKTIDX %lu) [%ld, %ld  <> %lu]",
                pkt_info.pktidx, wblk[0].block_num - 1, wblk[n_wblock-1].block_num + 1, pkt_blk_num);
          }
          else { // observing and first packet's timestamp is out of working range
            n_blks_rushed += pkt_blk_num - wblk[(n_wblock-1)/2].block_num;
            while(pkt_blk_num > wblk[(n_wblock-1)/2].block_num) { // only progress working range
              datablock_header = datablock_stats_header(&wblk[0]);
              hputu8(datablock_header, "PKTIDX", wblk[0].packet_idx);
              hputu8(datablock_header, "BLKSTART", wblk[0].packet_idx);
              hputu8(datablock_header, "BLKSTOP", wblk[1].packet_idx);
              // Finalize first working block
              finalize_block(wblk);
              // Update ndrop counter
              npacket_drop += wblk->ndrop;
              // hashpipe_info(thread_name, "Block dropped %d packets.", wblk->ndrop);
              // Shift working blocks
              block_stack_push(wblk, n_wblock);
              // Increment last working block
              increment_block(&wblk[n_wblock-1], wblk[n_wblock-1].block_num + 1);
              // Wait for new databuf data block to be free
              wait_for_block_free(&wblk[n_wblock-1], st, status_key);
            }
          }
      }
      observing = pkt_info.pktidx >= obs_start_seq_num && pkt_info.pktidx < obs_stop_seq_num;

      if (flag_reinit_blks) { // Reinitialise working blocks
        flag_reinit_blks = 0;
        // Re-init working blocks for block number of current packet's block,
        // and clear their data buffers
        for(wblk_idx=0; wblk_idx<n_wblock; wblk_idx++) {
          wblk[wblk_idx].pktidx_per_block = obs_info.pktidx_per_block;
          init_datablock_stats(wblk+wblk_idx, NULL, -1,
              wblk_idx,
              obs_info.pkt_per_block);
          wblk[wblk_idx].packet_idx = blk0_start_seq_num + wblk[wblk_idx].block_num * obs_info.pktidx_per_block;

          // also update the working blocks' headers
          datablock_header = datablock_stats_header(&wblk[wblk_idx]);
          hashpipe_status_lock_safe(st);
            memcpy(datablock_header, st->buf, HASHPIPE_STATUS_TOTAL_SIZE);
          hashpipe_status_unlock_safe(st);
        }
        hashpipe_info(thread_name, "Working block range now has PKTIDX range [%ld, %ld)", wblk[0].packet_idx, wblk[n_wblock-1].packet_idx + obs_info.pktidx_per_block);
      }

      // For each packet: process all packets
      #if ATA_IBV_FOR_PACKET_THREAD_COUNT > 1
        #pragma omp parallel for private (\
          p_pkt,\
          pkt_info,\
          p_payload,\
          blk0_relative_pkt_seq_num,\
          pkt_blk_num,\
          wblk_idx,\
          LATE_PKTIDX_flagged,\
          dest_feng_pktidx_offset\
        )\
        firstprivate (p_u8pkt, obs_info, PKT_OBS_FENG_flagged, PKT_OBS_SCHAN_flagged, PKT_OBS_NCHAN_flagged, PKT_OBS_PKTNTIME_flagged, PKT_OBS_PKTIDX_flagged)\
        reduction(min:obs_info_validity)\
        num_threads (ATA_IBV_FOR_PACKET_THREAD_COUNT)
        // The above `reduction` initialises each local variable as MAX>OBS_SEEMS_VALID, 
        //  and reduces (merges local variables) with min operator which preserves the
        //  negative values that indicate invalid obs_info
      #endif
      for(i=0; i < slots_per_block; i++) {
        if(PKT_OBS_FENG_flagged||PKT_OBS_SCHAN_flagged||PKT_OBS_NCHAN_flagged||PKT_OBS_PKTNTIME_flagged||PKT_OBS_PKTIDX_flagged){
          continue;
        }

        p_pkt = (struct ata_snap_ibv_pkt *)(p_u8pkt+i*slot_size);

        // Parse packet
        p_payload = (PKT_PAYLOAD_CP_T*)
          ata_snap_parse_ibv_packet(p_pkt, &pkt_info);

        // Get packet index and absolute block number for packet
        pkt_info.pktidx = pkt_info.pktidx - blk0_start_seq_num;

        // Only copy packet data and count packet if its wblk_idx is valid
        switch(check_pkt_observability_sans_idx(&obs_info, &pkt_info)){
          case PKT_OBS_OK:

            // Manage blocks based on pkt_blk_num
            // if (pkt_seq_num < blk0_start_seq_num && blk0_start_seq_num - pkt_seq_num < obs_info.pktidx_per_block){
            //   hashpipe_info(thread_name, "pkt_seq_num (%lu) < (%lu) blk0_start_seq_num", pkt_seq_num, blk0_start_seq_num);
            //   continue;
            // }
              
            // Once we get here, compute the index of the working block corresponding
            // to this packet.  The computed index may not correspond to a valid
            // working block!

            // Get packet's block number relative to the first block's starting index.
            pkt_blk_num = pkt_info.pktidx / obs_info.pktidx_per_block;
            wblk_idx = pkt_blk_num - wblk[0].block_num;

            if(0 <= wblk_idx && wblk_idx < n_wblock) {
              // Copy packet data to data buffer of working block
              dest_feng_pktidx_offset = ((PKT_PAYLOAD_CP_T*)
                datablock_stats_data(wblk+wblk_idx))
                + (pkt_info.pktidx%obs_info.pktidx_per_block)*time_byte_stride // offset for time// for TFP/_DP4A: 'Nants' is faster than 'Chans'
                + (pkt_info.pkt_schan-obs_info.schan)*channel_byte_stride
                + (pkt_info.feng_id*antenna_byte_stride); // offset for frequency
              
              #if ATA_IBV_TRANSPOSE_PACKET_THREAD_COUNT > 1
                #pragma omp parallel for\
                num_threads (ATA_IBV_TRANSPOSE_PACKET_THREAD_COUNT)
              #endif
              COPY_PACKET_PAYLOAD_FORLOOP(
                  dest_feng_pktidx_offset,
                  p_payload,
                  obs_info.pkt_nchan,
                  channel_byte_stride,
                  time_byte_stride);
              // Count packet for block and for processing stats
              thread_wblk_pkt_count[(omp_get_thread_num()*n_wblock) + wblk_idx] += 1;
            }
            else if(!LATE_PKTIDX_flagged){
              // Happens on the first packet that is outside of wblks' scope,
              // or if a packet arrives late, consider n_wblock++
              LATE_PKTIDX_flagged = 1;
              hashpipe_error(thread_name, "Packet ignored: determined wblk_idx = %d", wblk_idx);
            }
            break;
          case PKT_OBS_FENG:
            if(!PKT_OBS_FENG_flagged){
              obs_info_validity = OBS_INVALID_FENG;
              hashpipe_error(thread_name, 
                "Packet ignored: PKT_OBS_FENG\n\tfeng_id (%u) >= (%u) obs_info.nants",
                pkt_info.feng_id, obs_info.nants
              );
              PKT_OBS_FENG_flagged = 1;
            }
            break;
          case PKT_OBS_SCHAN:
            if(!PKT_OBS_SCHAN_flagged){
              obs_info_validity = OBS_INVALID_SCHAN;
              hashpipe_error(thread_name, 
                "Packet ignored: PKT_OBS_SCHAN\n\tpkt_schan (%d) < (%d) obs_info.schan",
                pkt_info.pkt_schan, obs_info.schan
              );
              PKT_OBS_SCHAN_flagged = 1;
            }
            break;
          case PKT_OBS_NCHAN:
            if(!PKT_OBS_NCHAN_flagged){
              obs_info_validity = OBS_INVALID_NCHAN;
              hashpipe_error(thread_name, 
                "Packet ignored: PKT_OBS_NCHAN\n\tpkt_chans [%d-%d] <> [%d-%d] obs_info.schan-obs_info.nchan",
                pkt_info.pkt_schan, pkt_info.pkt_schan + obs_info.pkt_nchan, obs_info.schan, obs_info.schan + obs_info.nchan
              );
              PKT_OBS_NCHAN_flagged = 1;
            }
            break;
          case PKT_OBS_PKTNTIME:
            if(!PKT_OBS_PKTNTIME_flagged){
              obs_info_validity = OBS_INVALID_PKTNTIME;
              hashpipe_error(thread_name, 
                "Packet ignored: PKT_OBS_PKTNTIME\n\tPKTNTIME %d != %d ATASNAP_DEFAULT_PKTNTIME",
                obs_info.pkt_ntime, ATASNAP_DEFAULT_PKTNTIME
              );
              PKT_OBS_PKTNTIME_flagged = 1;
            }
            break;
          case PKT_OBS_PKTIDX:
            if(!PKT_OBS_PKTIDX_flagged){
              obs_info_validity = OBS_INVALID_PKTIDX;
              hashpipe_error(thread_name, 
                "Packet ignored: PKT_OBS_PKTIDX\n\tPKTIDX-block0_pktidx %ld %% %d ATASNAP_DEFAULT_PKTNTIME != 0",
                pkt_info.pktidx, ATASNAP_DEFAULT_PKTNTIME
              );
              PKT_OBS_PKTIDX_flagged = 1;
            }
            break;
          default:
            obs_info_validity = OBS_INVALID;
            break;
        }
      } // end for each packet
      
      for(i=0; i < ATA_IBV_FOR_PACKET_THREAD_COUNT*n_wblock; i++){
        npacket += thread_wblk_pkt_count[i];
        wblk[i%n_wblock].npacket += thread_wblk_pkt_count[i];
      }
      memset(thread_wblk_pkt_count, 0, ATA_IBV_FOR_PACKET_THREAD_COUNT*n_wblock*sizeof(uint32_t));
    } // end if not idle
    else if (observing) {
      observing = 0;
    }

    // Mark input block free
    hpguppi_databuf_set_free(dbin, block_idx_in);

    // Handle 'min' reduced OBS_flags
    if(obs_info_validity == OBS_INVALID_FENG){
      hashpipe_status_lock_safe(st);
        hputs(st->buf, "OBSINFO", "INVALID FENG");
      hashpipe_status_unlock_safe(st);
    }
    else if(obs_info_validity == OBS_INVALID_SCHAN){
      hashpipe_status_lock_safe(st);
        hputs(st->buf, "OBSINFO", "INVALID SCHAN");
      hashpipe_status_unlock_safe(st);
    }
    else if(obs_info_validity == OBS_INVALID_NCHAN){
      hashpipe_status_lock_safe(st);
        hputs(st->buf, "OBSINFO", "INVALID NCHAN");
      hashpipe_status_unlock_safe(st);
    }
    else if(obs_info_validity == OBS_INVALID_PKTNTIME){
      hashpipe_status_lock_safe(st);
        hputs(st->buf, "OBSINFO", "INVALID PKTNCHAN");
      hashpipe_status_unlock_safe(st);
    }
    else if(obs_info_validity == OBS_INVALID_PKTIDX){
      hashpipe_status_lock_safe(st);
        hputs(st->buf, "OBSINFO", "INVALID PKTIDX !%%= ATASNAP_DEFAULT_PKTNTIME");
      hashpipe_status_unlock_safe(st);
    }
    else if(obs_info_validity == OBS_SEEMS_VALID) {
      obs_info_validity = OBS_VALID;
    }

    hashpipe_status_lock_safe(st);
    {
      hputu8(st->buf, "PKTIDX", pkt_info.pktidx);//wblk[0].packet_idx);

      // Update PKTIDX in status buffer if it is a new pkt_blk_num
      if(wblk[0].block_num != last_pkt_blk_num){
        last_pkt_blk_num = wblk[0].block_num;

          hputu8(st->buf, "BLKIDX", last_pkt_blk_num/obs_info.pktidx_per_block);
          hputu8(st->buf, "BLKSTART", wblk[0].packet_idx);
          hputu8(st->buf, "BLKSTOP", wblk[1].packet_idx);
      }
    }
    hashpipe_status_unlock_safe(st);
    
    if(wblk[n_wblock-1].npacket > 0) {
      // Time to advance the blocks!!!
      clock_gettime(CLOCK_MONOTONIC, &ts_stop_block);

      datablock_header = datablock_stats_header(&wblk[0]);
      hashpipe_status_lock_safe(st);
        memcpy(datablock_header, st->buf, HASHPIPE_STATUS_TOTAL_SIZE);
      hashpipe_status_unlock_safe(st);
      hputu8(datablock_header, "PKTIDX", wblk[0].packet_idx);
      hputu8(datablock_header, "BLKSTART", wblk[0].packet_idx);
      hputu8(datablock_header, "BLKSTOP", wblk[1].packet_idx);
      // Finalize first working block
      finalize_block(wblk);
      // Update ndrop counter
      npacket_drop += wblk->ndrop;
      // hashpipe_info(thread_name, "Block dropped %d packets.", wblk->ndrop);
      // Shift working blocks
      block_stack_push(wblk, n_wblock);
      // Increment last working block
      increment_block(&wblk[n_wblock-1], wblk[n_wblock-1].block_num + 1);
      // Wait for new databuf data block to be free
      wait_for_block_free(&wblk[n_wblock-1], st, status_key);

      blocks_per_second = 1e9/ELAPSED_NS(ts_start_block,ts_stop_block);
      clock_gettime(CLOCK_MONOTONIC, &ts_start_block);
    }

    // Update moving sum (for moving average)
    clock_gettime(CLOCK_MONOTONIC, &ts_free_input);
    fill_to_free_elapsed_ns = ELAPSED_NS(ts_stop_recv, ts_free_input);
    // Add new value, subtract old value
    fill_to_free_moving_sum_ns +=
        fill_to_free_elapsed_ns - fill_to_free_block_ns[block_idx_in];
    // Store new value
    fill_to_free_block_ns[block_idx_in] = fill_to_free_elapsed_ns;

#if 0
    fprintf(stderr, "blkin %d fill at %ld free +%ld ns (%d packets)\n",
        block_idx_in,
        ELAPSED_NS(ts_stop_recv, ts_free_input), njobs);
#endif

    // Advance to next input block
    block_idx_in = (block_idx_in + 1) % dbin->header.n_block;

    // Will exit if thread has been cancelled
    pthread_testcancel();
  } // end main loop

  hashpipe_info(thread_name, "exiting!");
  pthread_exit(NULL);

  return NULL;
}

static hashpipe_thread_desc_t thread_desc = {
    name: "hpguppi_ata_ibv_payload_order_thread",
    skey: "NETSTAT",
    init: init,
    run:  run,
    ibuf_desc: {hpguppi_input_databuf_create},
    obuf_desc: {ATA_IBV_OUT_DATABUF_CREATE}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(&thread_desc);
}

// vi: set ts=2 sw=2 et :
