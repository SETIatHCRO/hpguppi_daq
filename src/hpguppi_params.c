/* hpguppi_params.c
 *
 * Routines to read/write basic system parameters.
 * Use PSRFITS style keywords as much as possible.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <hashpipe.h>

#include "hpguppi_params.h"
#include "hpguppi_time.h"

#include "fitshead.h"
#include "psrfits.h"
#include "slalib.h"

#ifndef DEGTORAD
#define DEGTORAD 0.017453292519943295769236907684886127134428718885417
#endif
#ifndef RADTODEG
#define RADTODEG 57.29577951308232087679815481410517033240547246656
#endif
#ifndef SOL
#define SOL 299792458.0
#endif

#define DEBUGOUT 0

#define get_dbl(key, param, def) {                                      \
        if (hgetr8(buf, (key), &(param))==0) {                          \
            if (DEBUGOUT)                                               \
                printf("Warning:  %s not in status shm!\n", (key));     \
            (param) = (def);                                            \
        }                                                               \
    }

#define get_int(key, param, def) {                                      \
        if (hgeti4(buf, (key), &(param))==0) {                          \
            if (DEBUGOUT)                                               \
                printf("Warning:  %s not in status shm!\n", (key));     \
            (param) = (def);                                            \
        }                                                               \
    }

#define get_uint(key, param, def) {                                      \
        if (hgetu4(buf, (key), &(param))==0) {                          \
            if (DEBUGOUT)                                               \
                printf("Warning:  %s not in status shm!\n", (key));     \
            (param) = (def);                                            \
        }                                                               \
    }

#define get_lon(key, param, def) {                                      \
        {                                                               \
            double dtmp;                                                \
            if (hgetr8(buf, (key), &dtmp)==0) {                         \
                if (DEBUGOUT)                                           \
                    printf("Warning:  %s not in status shm!\n", (key)); \
                (param) = (def);                                        \
            } else {                                                    \
                (param) = (long long)(rint(dtmp));                      \
            }                                                           \
        }                                                               \
    }

#define get_str(key, param, len, def) {                                 \
        if (hgets(buf, (key), (len), (param))==0) {                     \
            if (DEBUGOUT)                                               \
                printf("Warning:  %s not in status shm!\n", (key));     \
            strcpy((param), (def));                                     \
        }                                                               \
    }

#define exit_on_missing(key, param, val) {                              \
        if ((param)==(val)) {                                           \
            char errmsg[100];                                           \
            sprintf(errmsg, "%s is required!\n", (key));                \
            hashpipe_error(__FUNCTION__, errmsg);                       \
            exit(1);                                                    \
        }                                                               \
    }


#if 0
// Return the beam FWHM in degrees for obs_freq in MHz
// and dish_diam in m
double beam_FWHM(double obs_freq, double dish_diam)
{
    double lambda = SOL/(obs_freq*1e6);
    return 1.2 * lambda / dish_diam * RADTODEG;
}
#endif // 0


// Any GB-specific derived parameters go here
static void set_obs_params_gb(char *buf,
                              struct hpguppi_params *g,
                              struct psrfits *p) {

    // TODO could double-check telescope name first
    // Load the values from status memory or fill-in with default values
    get_int("FD_HAND", p->hdr.fd_hand, -1);
    get_dbl("FD_SANG", p->hdr.fd_sang, 45.0);
    get_dbl("FD_XYPH", p->hdr.fd_xyph, 0.0);
    get_int("BE_PHASE", p->hdr.be_phase, -1);
    get_dbl("BEAMFWHM", p->hdr.beam_FWHM, 65.0);
}

// Read networking parameters
void hpguppi_read_net_params(char *buf, struct hashpipe_udp_params *u) {
    get_str("DATAHOST", u->sender, sizeof(u->sender), "bee2-10");
    get_int("DATAPORT", u->port, 50000);
    get_str("PKTFMT", u->packet_format, sizeof(u->packet_format), "GUPPI");
    if (strncmp(u->packet_format, "PARKES", 6)==0)
        u->packet_size = 2056;
    else if (strncmp(u->packet_format, "1SFA", 4)==0)
        u->packet_size = 8224;
    else if (strncmp(u->packet_format, "FAST4K", 6)==0)
        u->packet_size = 4128;
    else if (strncmp(u->packet_format, "SHORT", 5)==0)
        u->packet_size = 544;
    else
        u->packet_size = 8208;
}

// Read networking parameters for packet sockets.  Same as for UDP sockets,
// though BINDHOST rather than DATAHOST must be used and it should be a local
// interface name (e.g. eth4).  Also BINDPORT rather than DATAPORT.
void hpguppi_read_pktsock_params(char *buf, struct hpguppi_pktsock_params *p)
{
    get_str("BINDHOST", p->ifname, 80, "eth4");
    get_int("BINDPORT", p->port, 60000);
    get_int("OBSSCHAN", p->obsschan, 0);
    get_int("OBSNCHAN", p->obsnchan, 64);
    get_str("PKTFMT", p->packet_format, 32, "1SFA");
    if (strncmp(p->packet_format, "PARKES", 6)==0)
        p->packet_size = 2056;
    else if (strncmp(p->packet_format, "1SFA", 4)==0)
        p->packet_size = 8224;
    else if (strncmp(p->packet_format, "FAST4K", 6)==0)
        p->packet_size = 4128;
    else if (strncmp(p->packet_format, "SHORT", 5)==0)
        p->packet_size = 544;
    else if (strncmp(p->packet_format, "ATASNAPV", 8)==0){
        int bitspersample;
        get_int("NBITS", bitspersample, 4);
        //               payload : complex4[nchan, 16, 2]      +  header 
        p->packet_size = p->obsnchan*16*2*(bitspersample*2/8) + 16;
    }
    else
        p->packet_size = 8208;
}

#if 0
/* Some code just needs a simple way to get the obs mode string */
void guppi_read_obs_mode(const char *buf, char *mode) {
    get_str("OBS_MODE", mode, 24, "Unknown");
}
#endif // 0

// Read a status buffer all of the key observation paramters
// Expectedly called after hpguppi_read_obs_params
void hpguppi_read_subint_params(char *buf,
                                struct hpguppi_params *g,
                                struct psrfits *p)
{
    // Parse packet size, # of packets, etc.
    get_lon("PKTIDX", g->packetindex, -1L);
    get_int("NPKT", g->n_packets, 0);
    get_int("NDROP", g->n_dropped, 0);
    get_dbl("DROPAVG", g->drop_frac_avg, 0.0);
    get_dbl("DROPTOT", g->drop_frac_tot, 0.0);
    if (g->packetsize==0)
        get_int("PKTSIZE", g->packetsize, 0);

    if (g->n_packets>0)
        g->drop_frac = (double) g->n_dropped / (double) g->n_packets;
    else
        g->drop_frac = 0.0;

    // Valid obs start time
    get_int("STTVALID", g->stt_valid, 0);

    // Observation params
    get_dbl("AZ", p->sub.tel_az, 0.0);
    if (p->sub.tel_az < 0.0) p->sub.tel_az += 360.0;
    get_dbl("ZA", p->sub.tel_zen, 0.0);
    get_dbl("RA", p->sub.ra, 0.0);
    get_dbl("DEC", p->sub.dec, 0.0);

    // Backend HW parameters
    get_int("ACC_LEN", g->decimation_factor, 0);
    get_int("NBITSADC", g->n_bits_adc, 8);
    get_int("PFB_OVER", g->pfb_overlap, 4);
    get_int("CODD", g->coherent, 0);

    // Check fold mode
    int fold=0;
    if (strstr("PSR", p->hdr.obs_mode)!=NULL) { fold=1; }
    if (strstr("CAL", p->hdr.obs_mode)!=NULL) { fold=1; }

    // Fold-specifc stuff
    if (fold) {
        get_dbl("TSUBINT", p->sub.tsubint, 0.0);
        get_dbl("OFFS_SUB", p->sub.offs, 0.0);
        get_int("NPOLYCO", p->fold.n_polyco_sets, 0);
    } else {
        int bytes_per_dt = p->hdr.nchan * p->hdr.npol * p->hdr.nbits / 8;
        p->sub.offs = p->hdr.dt *
            (double)(g->packetindex * g->packetsize / bytes_per_dt)
            + 0.5 * p->sub.tsubint;
        p->fold.n_polyco_sets = 0;
    }

    { // MJD and LST calcs
        int imjd, smjd, lst_secs;
        double offs, mjd;
        get_current_mjd(&imjd, &smjd, &offs);
        mjd = (double) imjd + ((double) smjd + offs) / 86400.0;
        get_current_lst(mjd, &lst_secs);
        p->sub.lst = (double) lst_secs;
    }

    // Until we need them...
    p->sub.feed_ang = 0.0;
    p->sub.pos_ang = 0.0;
    p->sub.par_ang = 0.0;

    // Galactic coords
    slaEqgal(p->sub.ra*DEGTORAD, p->sub.dec*DEGTORAD,
             &p->sub.glon, &p->sub.glat);
    p->sub.glon *= RADTODEG;
    p->sub.glat *= RADTODEG;
}


// Read a status buffer all of the key observation paramters
// Call once at observation start, then call lighter-weight 
// hpguppi_read_subint_params
void hpguppi_read_obs_params(char *buf,
                             struct hpguppi_params *g,
                             struct psrfits *p)
{
    char base[200], dir[200], banknam[64];

    // Software data-stream modification params
    get_int("DS_TIME", p->hdr.ds_time_fact, 1); // Time down-sampling
    get_int("DS_FREQ", p->hdr.ds_freq_fact, 1); // Freq down-sampling
    get_int("ONLY_I", p->hdr.onlyI, 0);         // Only output Stokes I

    // PIPERBLK
    get_int("PIPERBLK", g->packets_per_block, 0);
    if (g->packets_per_block==0)
    {
        get_int("BLOCSIZE", g->packets_per_block, 0);
        get_int("PKTSIZE", g->packetsize, 0);
        if (g->packetsize>0)
            g->packets_per_block /= g->packetsize;
    }

    // Freq, BW, etc.
    get_dbl("OBSFREQ", p->hdr.fctr, 0.0);
    get_dbl("OBSBW", p->hdr.BW, 0.0);
    exit_on_missing("OBSBW", p->hdr.BW, 0.0);
    get_int("OBSNCHAN", p->hdr.nchan, 0);
    exit_on_missing("OBSNCHAN", p->hdr.nchan, 0);
    get_int("NPOL", p->hdr.npol, 0);
    exit_on_missing("NPOL", p->hdr.npol, 0);
    get_int("NBITS", p->hdr.nbits, 0);
    exit_on_missing("NBITS", p->hdr.nbits, 0);
    get_dbl("TBIN", p->hdr.dt, 0.0);
    exit_on_missing("TBIN", p->hdr.dt, 0.0);
    get_dbl("CHAN_BW", p->hdr.df, 0.0);
    if (p->hdr.df==0.0) p->hdr.df = p->hdr.BW/p->hdr.nchan;
    get_dbl("SCANLEN", p->hdr.scanlen, 0.0);
    get_int("NRCVR", p->hdr.rcvr_polns, 2);
    p->hdr.orig_df = p->hdr.df;
    p->hdr.orig_nchan = p->hdr.nchan;

    // Observation information
    get_str("TELESCOP", p->hdr.telescope, 24, "GBT");
    get_str("OBSERVER", p->hdr.observer, 24, "Unknown");
    get_str("SRC_NAME", p->hdr.source, 24, "Unknown");
    get_str("FRONTEND", p->hdr.frontend, 24, "Unknown");
    get_str("BACKEND", p->hdr.backend, 24, "GUPPI");
    get_str("PROJID", p->hdr.project_id, 24, "Unknown");
    get_str("FD_POLN", p->hdr.poln_type, 8, "Unknown");
    get_str("POL_TYPE", p->hdr.poln_order, 16, "Unknown");
    get_int("SCANNUM", p->hdr.scan_number, 1);
    get_str("DATADIR", dir, 200, ".");
    get_str("BANKNAM", banknam, sizeof(banknam), ".");
    if (strcmp(p->hdr.poln_order, "AA+BB")==0 ||
        strcmp(p->hdr.poln_order, "INTEN")==0)
        p->hdr.summed_polns = 1;
    else
        p->hdr.summed_polns = 0;
    get_str("TRK_MODE", p->hdr.track_mode, 16, "Unknown");
    get_str("RA_STR", p->hdr.ra_str, 16, "Unknown");
    get_str("DEC_STR", p->hdr.dec_str, 16, "Unknown");
    // Should set other cal values if CAL_MODE is on
    get_str("CAL_MODE", p->hdr.cal_mode, 8, "Unknown");
    if (!(strcmp(p->hdr.cal_mode, "OFF")==0)) {  // Cals not off
        get_dbl("CAL_FREQ", p->hdr.cal_freq, 25.0);
        get_dbl("CAL_DCYC", p->hdr.cal_dcyc, 0.5);
        get_dbl("CAL_PHS", p->hdr.cal_phs, 0.0);
    }
    get_str("OBS_MODE", p->hdr.obs_mode, 24, "Unknown");

    // Fold mode specific stuff
    int fold=0;
    get_int("NBIN", p->fold.nbin, 256);
    get_dbl("TFOLD", p->fold.tfold, 30.0);
    get_str("PARFILE", p->fold.parfile, 256, "");
    if (strstr("FOLD", p->hdr.obs_mode)!=NULL) { fold=1; }
    if (strstr("PSR", p->hdr.obs_mode)!=NULL) { fold=1; }
    if (strstr("CAL", p->hdr.obs_mode)!=NULL) { fold=1; }
    if (fold)
        p->hdr.nbin = p->fold.nbin;
    else
        p->hdr.nbin = 1;

    // Coherent dedispersion params
    get_int("CODD", g->coherent, 0);
    get_lon("OBSSTART", g->start_pkt, -1);
    if(g->start_pkt == -1){
        get_lon("PKTSTART", g->start_pkt, 0);
    }
    get_int("OVERLAP", p->dedisp.overlap, 0);
    get_dbl("CHAN_DM", p->hdr.chan_dm, 0.0);

    { // Start time, MJD
        int mjd_d, mjd_s;
        double mjd_fs;
        get_int("STT_IMJD", mjd_d, 0);
        get_int("STT_SMJD", mjd_s, 0);
        get_dbl("STT_OFFS", mjd_fs, 0.0);
        p->hdr.MJD_epoch = (long double) mjd_d;
        p->hdr.MJD_epoch += ((long double) mjd_s + mjd_fs) / 86400.0;
        p->hdr.start_day = mjd_d;
        p->hdr.start_sec = mjd_s + mjd_fs;
    }

    // Set the base filename
    int i;
    char backend[24];
    strncpy(backend, p->hdr.backend, 24);
    for (i=0; i<24; i++) {
        if (backend[i]=='\0') break;
        backend[i] = tolower(backend[i]);
    }
    if (strstr("CAL", p->hdr.obs_mode)!=NULL) {
        sprintf(base, "%s_%05d_%s_%04d_cal", backend, p->hdr.start_day,
                p->hdr.source, p->hdr.scan_number);
    } else {
        // base is BACKEND_MJD_SEC_BLK_SRC_SCAN.
        sprintf(base, "%s_%05d_%05d_%06lld_%s_%04d", backend,
                p->hdr.start_day, (int)(p->hdr.start_sec), g->start_pkt/g->packets_per_block,
                p->hdr.source, p->hdr.scan_number);
    }
#ifdef NO_PROJECT_DIR
    sprintf(p->basefilename, "%s/%s", dir, base);
#else
    // Use a $DATADIR/$PROJID/$BACKEND/$BANK prefix for files
    if (strnlen(banknam, sizeof(banknam)) < 1)
        snprintf(banknam, sizeof(banknam), ".");
    sprintf(p->basefilename, "%.72s/%s/%s%s/%c/%.71s", dir, p->hdr.project_id,
            p->hdr.backend, g->coherent ? "_CODD" : "",
            banknam[strnlen(banknam, sizeof(banknam))-1], base);
#endif
    { // Date and time of start
        int YYYY, MM, DD, h, m, s_integ_int, s_frac_int;
        double s, s_integ, s_frac;
        datetime_from_mjd(p->hdr.MJD_epoch, &YYYY, &MM, &DD, &h, &m, &s);
        s_frac = modf(s, &s_integ);
        s_integ_int = abs((int)s_integ);
        s_frac_int = abs((int)(s_frac*1000));
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wformat-overflow="
        sprintf(p->hdr.date_obs, "%04d-%02d-%02dT%02d:%02d:%02d.%03d",
                YYYY % 10000, MM % 100, DD % 100, h % 24, m % 60, s_integ_int % 100, s_frac_int % 1000);
        #pragma GCC diagnostic pop
    }

    // TODO: call telescope-specific settings here
    // Eventually make this depend on telescope name
    set_obs_params_gb(buf, g, p);

    // Now bookkeeping information
    {
        int ii, jj, kk;
        int bytes_per_dt = p->hdr.nchan * p->hdr.npol * p->hdr.nbits / 8;
        char key[17];
        double offset, scale, dtmp;
        long long max_bytes_per_file;

        get_int("BLOCSIZE", p->sub.bytes_per_subint, 0);
        p->hdr.nsblk = p->sub.bytes_per_subint / bytes_per_dt;
        p->sub.FITS_typecode = TBYTE;
        p->sub.tsubint = p->hdr.nsblk * p->hdr.dt;
        if (fold) {
            //p->hdr.nsblk = 1;
            p->sub.FITS_typecode = TFLOAT;
            get_dbl("TSUBINT", p->sub.tsubint, 0.0);
            p->sub.bytes_per_subint = sizeof(float) * p->hdr.nbin *
                p->hdr.nchan * p->hdr.npol;
            //max_bytes_per_file = PSRFITS_MAXFILELEN_FOLD * 1073741824L;
            max_bytes_per_file = PSRFITS_MAXFILELEN_FOLD * 1073741824L / 8L;
        } else {
            max_bytes_per_file = PSRFITS_MAXFILELEN_SEARCH * 1073741824L;
        }
        // Will probably want to tweak this so that it is a nice round number
        if (p->sub.bytes_per_subint!=0)
            p->rows_per_file = p->hdr.ds_freq_fact * p->hdr.ds_time_fact *
                (p->hdr.onlyI ? 4 : 1) * max_bytes_per_file /
                p->sub.bytes_per_subint;

        // Free the old ones in case we've changed the params
        hpguppi_free_psrfits(p);

        // Allocate the subband arrays
        p->sub.dat_freqs = (float *)malloc(sizeof(float) * p->hdr.nchan);
        p->sub.dat_weights = (float *)malloc(sizeof(float) * p->hdr.nchan);
        // The following correctly accounts for the middle-of-bin FFT offset
        dtmp = p->hdr.fctr - 0.5 * p->hdr.BW;
        double foffs;
        int is_codd;
        get_int("CODD", is_codd, 0);
        if (is_codd) { foffs=0.5; }
        else { foffs=0.0; }
        for (ii = 0 ; ii < p->hdr.nchan ; ii++) {
            p->sub.dat_freqs[ii] = dtmp + (ii+foffs) * p->hdr.df;
            p->sub.dat_weights[ii] = 1.0;
        }
        // Explicitly weight the DC and Nyquist channels zero
        // because of how power is split between them
        // XXX this needs to be changed for coherent dedisp..
        //p->sub.dat_weights[0] = 0.0;

        p->sub.dat_offsets = (float *)malloc(sizeof(float) *
                                             p->hdr.nchan * p->hdr.npol);
        p->sub.dat_scales = (float *)malloc(sizeof(float) *
                                            p->hdr.nchan * p->hdr.npol);
        for (ii = 0 ; ii < p->hdr.npol ; ii++) {
            sprintf(key, "OFFSET%d", ii);
            get_dbl(key, offset, 0.0);
            sprintf(key, "SCALE%d", ii);
            get_dbl(key, scale, 1.0);
            for (jj = 0, kk = ii*p->hdr.nchan ; jj < p->hdr.nchan ; jj++, kk++) {
                p->sub.dat_offsets[kk] = offset;
                p->sub.dat_scales[kk] = scale;
            }
        }
    }

    // Read information that is appropriate for the subints
    hpguppi_read_subint_params(buf, g, p);
    p->hdr.azimuth = p->sub.tel_az;
    p->hdr.zenith_ang = p->sub.tel_zen;
    p->hdr.ra2000 = p->sub.ra;
    p->hdr.dec2000 = p->sub.dec;
    p->hdr.start_lst = p->sub.lst;
    p->hdr.feed_angle = p->sub.feed_ang;
}

void hpguppi_free_psrfits(struct psrfits *p)
{
    // Free any memory allocated in to the psrfits struct
    if (p->sub.dat_freqs) free(p->sub.dat_freqs);
    if (p->sub.dat_weights) free(p->sub.dat_weights);
    if (p->sub.dat_offsets) free(p->sub.dat_offsets);
    if (p->sub.dat_scales) free(p->sub.dat_scales);
}

// Read direct I/O mode.  This comes from the DIRECTIO status buffer key.  It
// is interpreted as a number.  Defined and non-zero numeric (e.g. 1) means to
// use Direct I/O.  Undefined or zero or non-numeric means do NOT use Direct
// I/O.
int hpguppi_read_directio_mode(char *buf)
{
    int directio = 0;
    get_int("DIRECTIO", directio, 0);
    return directio;
}

// Read PIPERBLK.  PIPERBLK is the amount PKTIDX increments per block.
// If zero or missing, it must be inferred from the change in PKTIDX from one
// block to the next.
unsigned int hpguppi_read_piperblk(char *buf)
{
    int piperblk = 0;
    get_int("PIPERBLK", piperblk, 0);
    return piperblk;
}

// Calculate the largest power of two number of time samples that fit in
// max_block_size block size for a given number of channels.
int calc_ntime_per_block(int max_block_size, int obsnchan)
{
    int ntime_per_block = max_block_size / (4 * obsnchan);

    // Set ntime_per_block to closest power of 2 less than or equal to
    // ntime_per_block.
    unsigned int i = 0;
    while(ntime_per_block != 1) {
        ntime_per_block >>= 1;
        i++;
    }
    ntime_per_block <<= i;

    return ntime_per_block;
}
