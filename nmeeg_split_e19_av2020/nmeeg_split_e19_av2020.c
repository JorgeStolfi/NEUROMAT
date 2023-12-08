#define PROG_NAME "nmeeg_split_e19_av2020"
#define PROG_DESC "Parses and cleans the EEG session data file of A. Valencio's 2020-02 project; splits into single-probe files."
#define PROG_VERS "2020-02-27"

#define nmeeg_split_e19_av2020_C_COPYRIGHT \
  "Copyright © 2020 by the State University of Campinas (UNICAMP)"

/* Last edited on 2023-11-02 06:29:52 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    [ -firstBlockID {FBLOCK} ] \\\n" \
  "    [ -firstProbeID {FPROBE} ] \\\n" \
  "    -padFrames {PROBE_PAD} \\\n" \
  "    -outDir {OUT_DIR} \\\n" \
  "    " argparser_help_info_HELP " \\\n" \
  "    < {INFILE}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program reads from standard input a 19-electrode EEG data file" \
  " from the Deep Brain Stimulation (DBS) experiment by Arthur Valencio et al. (collected in 2020-02), filters out certain noise components, and" \
  " extracts from it the segments that correspond to individual experimental" \
  " probes, as separate files.\n" \
  "\n" \
  "INPUT FORMAT\n" \
  "  The input file should be a full-session EEG recording, sampled at 500 Hz, converted to \".txt\" format.\n" \
  "\n" \
  "  The program assumes that the input file contains a certain" \
  " number {NT} of data frames, each with the same number {NC = 21} of" \
  " channels; including {NE = NC-1 = 20} electrode potential channels, followed" \
  " by one screen brightness sensor (SBS) channel.  The electrode channels include 19 actual" \
  " electrodes and a synthetic one, for the purpose of monitoring" \
  " the filtering and extraction operations.  The SBS channel" \
  " identifies different segments of the recording.\n" \
  "\n" \
  "  The program assumes that the EEG file has the following structure:\n" \
  "\n" \
  "    The whole file is a SESSION.\n" \
  "\n" \
  "    The SESSION consists of two PHASES, 'DBS_off' and 'DBS_on', interleaved" \
  " with three SETUPS.\n" \
  "\n" \
  "    During each SETUP, the SBS channel has value 0.5.\n" \
  "\n" \
  "    Each PHASE consists of three BLOCKS interleaved with two REST periods.\n" \
  "\n" \
  "    Each REST period lasts about 30 seconds and has the SBS channel" \
  " at about 0.5.\n" \
  "\n" \
  "    Each BLOCK consists of 10 PROBES interleaved with 9 GAPS.\n" \
  "\n" \
  "    Each GAP lasts about 3 seconds, and, like a REST, has the SBS channel" \
  " at about 0.5.\n" \
  "\n" \
  "    Each PROBE lasts about 30 seconds and can be of two types: 'human static', with" \
  " the SBS channel at 0, or 'human moving', with SBS channel at 1.\n" \
  "\n" \
  "  The program assumes that the original SBS channel (the voltage out of the" \
  " device, from 0 to about -0.26 volts) has been cleaned somehow and" \
  " has been mapped to the three exact values 0.0. 0.5, and 1.0.  In" \
  " particular, the blinking used to identify the SETUP segments has been turned to steady 0.5.\n" \
  "\n" \
  "  The program will first extract from the file two sub-files, each containing" \
  " the frames for one PHASE of the SESSION, plus 2 times {PROBE_PAD} frames from the" \
  " GAP, REST, or SETUP intervals on either side.  The program will then filter" \
  " each electrode signal in each phase, with a set of notch filters that remove" \
  " certain peak frequencies, and a band-pass filter.\n" \
  "\n" \
  "  The program will then split the filtered PHASE sub-files into a separate file for each PROBE" \
  " segment.  The PROBE file will have exactly {PROBE_PAD} frames from the" \
  " previous GAP, REST, or SETUP intervals, then the PROBE segment, then enough frames" \
  " from the following GAP, REST, or SETUP to complete exactly {PROBE_PAD + PROBE_NFR + PROBE_PAD} frames.\n" \
  "\n" \
  "OUTPUT FILES\n" \
  "  The output files are written with" \
  " names \"{OUT_DIR}/s{SSS}_dbs{D}.txt\" and \"{OUT_DIR}/s{SSS}_dbs{D}_r{BBB}{RR}.txt\" where {SSS}" \
  " is the three-digit subject ID number (obtained from the input file's header), {D} is" \
  " the PHASE index (0 or 1 according to whether the DBS is turned off or on, respectively), {BBB} is" \
  " a three-digit sequential block number within the file, starting with {FBLOCK}, and {RR} is" \
  " a two-digit sequential probe number within the block, starting with {FPROBE}.  (See" \
  " the options \"-firstBlockID\" and \"-firstProbeID\" below.)\n" \
  "\n" \
  "  Each output frame contains the same number {NC = NE+1} channels, comprising the same {NE} electrode" \
  " channels as in the input followed by the SBS channel.\n" \
  "\n" \
  "  Each output file begins with header records that specify some relevant" \
  " parameters such as number of frames and channels, sampling frequency, etc..\n" \
  "\n" \
  "OPTIONS\n" \
  "  -firstBlockID {FBLOCK}\n" \
  "    This optional argument specifies the ID number of the first block" \
  " of probes in the input file, among all probes of the subject.  If" \
  " omitted, the program assumes {FBLOCK=1}.\n" \
  "\n" \
  "  -firstProbeID {FPROBE}\n" \
  "    This optional argument specifies the ID number of the first probe\n" \
  " in each block of probes.  If omitted, assumes {FPROBE=1}.\n" \
  "\n" \
  "  -padFrames {PROBE_PAD}\n" \
  "    This mandatory argument specifies the number of non-PROBE frames (GAP, REST, or SETUP) to extract" \
  " before and after each PROBE segment.\n" \
  "\n" \
  "  -outDir {OUT_DIR}\n" \
  "    This mandatory argument specifies the directory where the" \
  " output files will be placed (see the OUTPUT FILES section above). The" \
  " directory must exist and must be writable.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  neuromat_eeg_plot_signals.sh(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2020-02-27 by Jorge Stolfi, IC-UNICAMP, based on {nmeeg_split_e128_gh2013.c}.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2020-02-27 Created (J.Stolfi).\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " nmeeg_split_e19_av2020_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <math.h>

#include <fftw3.h>

#include <argparser.h>
#include <vec.h>
#include <affirm.h>
#include <jsfile.h>
#include <jsmath.h>
#include <jsstring.h>

#include <neuromat_eeg.h>
#include <neuromat_filter.h>
#include <neuromat_eeg_io.h>
#include <neuromat_eeg_header.h>
#include <neuromat_eeg_source.h>

typedef struct nes_options_t
  { 
    int firstBlockID;   /* External number of first probe block found in input file. */
    int firstProbeID;   /* External number of first probe in each block. */
    int padFrames;    /* Desired number of padding data frames around each extracted PROBE. */
    char *outDir;     /* Directory for all output files. */
  } nes_options_t;
  /* Arguments from command line. */
  
#define nes_SAMPLING_FREQ 500.0
  /* Sampling frequency (Hz; samples per second). */
  
#define nes_NUM_INPUT_CHANNELS 21
  /* Number of channels.  Includes the real EEG electrodes, a sythetic signal 
    to monitor the effects of filtering, and the (cleaned) SBS signal. */

#define nes_NUM_SCALP_ELECTRODES 19
  /* Number of real EEG electrodes. */

#define nes_SUBJECT_ID 2
  /* ID number of the only subject. */
  
#define nes_NUM_INPUT_FRAMES 1466500
  /* Expected number of frames in whole SESSION file. */
  
#define nes_NUM_PHASES_PER_SESSION 2
  /* Expected number of PHASES in SESSION file. */
  
#define nes_NUM_BLOCKS_PER_PHASE 3
  /* Expected number of BLOCKS per PHASE in file. */

#define nes_NUM_PROBES_PER_BLOCK 10
  /* Expected number of PROBES in each BLOCK in file. */

#define nes_NUM_FRAMES_PER_PROBE (30*500)
  /* Expected number of data (stimulus) frames per PROBES. */

int main(int argc, char **argv);
  /* Main prog. */
  
nes_options_t *nes_parse_options(int argc, char **argv);
  /* Parses the command line arguments. */
  
void nes_get_PROBES
  ( int32_t nt, 
    int32_t nc, 
    double **smp, 
    double fsmp,
    int32_t ic_sbs,
    int32_t np, 
    int32_t nb, 
    int32_t nr,
    int32_t kt_ini[], /* (OUT) */
    int32_t kt_fin[]  /* (OUT) */
  );
  /* Identifies the initial and final frame indices of each PROBE in the
    EEG dataset {smp[0..nt-1][0..nc-1]}.
    
    The structure of the dataset is assumed to be marked by the SBS
    channel (column {ic_sbs} of the array). Each PHASE is expectd to
    have {nb} BLOCKS of {nr} PROBES, during which the SBS signal is 0.0
    or 1.0, surrounded by SETUP, REST, or GAP intervals, in which the
    SBS signal is assumed to be 0.5.
    
    The indices of the first and last frame of PROBE with index {ir} of
    BLOCK with index {ib} and PHASE with index {ip} is returned
    {kt_ini[j]} and ,kt_fin[j]}, respectively, where {j} is {ir + (ib +
    ip*nb)*nr}. */

void nes_PHASE_filter
  ( int32_t nt, 
    int32_t nc, 
    double **smp, 
    double fsmp,
    int32_t ic_sbs,
    int32_t ip,
    int32_t it_ini,
    int32_t it_fin
  );
  /* Filters each data channel of PHASE with index {ip}.
    Assumes that PHASE frames (including padding) are
    rows {it_ini..it_fin} of {smp}, and the data channels (including
    the synthetic one) are columns {0..nc-1}, except the SBS signal
    in column {ic_sbs}. */

void nes_signal_filter(int32_t ns, double fsmp, double x[], double X[], double F[], fftw_plan pdir, fftw_plan pinv);
  /* Filters the signal {x[0..nt-1]} with the Fourier/Hartley weights
    {F[0..ns-1]}, assuming that the sampling frequency is {fsmp} The
    FFTW plans {pdir,pinv} must map {x[0..ns-1]} to/from {X[0..ns-1]}.
    Better allocate {x} and {X} with extra elems to be sure. */
    
double *nes_make_filter(int32_t nt, double fsmp);
  /* Creates a real symmetric (phase-preserving) Fourier filter for {nt} data samples,
   assuming the sampling frequency {fsmp}.  The filter removes certain narrow
   frequency bands due to DBS interfernce, and then applies a bandpass filter
   with gradual cutoffs at 1 Hz and 50 Hz. */

void nes_PHASE_write
  ( char *outDir,
    int32_t nt, 
    int32_t nc, 
    double **smp, 
    double fsmp,
    int32_t ic_syn,
    int32_t ic_sbs,
    int32_t subject,
    int32_t ip,
    int32_t it_ini,
    int32_t it_fin
  );
  /* Writes the EEG data {smp[it_ini..it_fin][0..nc-1]} to file
    "{outDir}/s{SSS}_p{P}.txt", where {SSS} is the {subject} (3 digits), {P}
    is the PHASE index {ip} (0 or 1, 1 digit), all zero-padded. All {nc} channels
    are written, including the synthetic signal (column {ic_syn}) and 
    the SBS signal (column {ic_sbs}). */
  
void nes_PROBE_write
  ( char *outDir,
    int32_t nt, 
    int32_t nc, 
    double **smp, 
    double fsmp,
    int32_t ic_syn,
    int32_t ic_sbs,
    int32_t subject,
    int32_t ip,
    int32_t ib,
    int32_t ir,
    int32_t it_ini,
    int32_t it_fin
  );
  /* Writes the EEG data {smp[it_ini..it_fin][0..nc-1]} to file
    "{outDir}/s{SSS}_p{P}_b{B}_r{RR}.txt", where {SSS} is the {subject} ID (3 digits), {P}
    is the PHASE index {ip} (0 or 1, 1 digit), {B} is the BLOCK ID ({ib+1}, 1 digit),
    {RR} is the PROBE ID ({ir+1}, 2 digits), all zero-padded. EXCLUDES the
    synthetic signal (column {ic_syn}) but writes the SBS signal
    (column {ic_sbs}). */
  
void nes_header_write
  ( FILE *wr, 
    int32_t nt, 
    int32_t nc, 
    double fsmp,
    int32_t it_step,
    int32_t ic_syn,
    int32_t ic_sbs,
    int32_t subject,
    char *type,
    int32_t it_ini,
    int32_t it_fin
  );
  /* Writes a file header in the format described in {neuromat_eeg_header.h}.
    The {type} should be "DBS0" or "DBS1" for the two phases. */
  

int main(int argc, char **argv)
  {
    nes_options_t *o = nes_parse_options(argc, argv);
    /* bool_t debug = TRUE; */
    /* bool_t debugFrame = FALSE; */ /* Print debug for each frame read. */
    
    /* Print some options: */
    fprintf(stderr, "first experimental block in input file is number %03d\n", o->firstBlockID);
    fprintf(stderr, "first probe of each block is number %02d\n", o->firstProbeID);
    fprintf(stderr, "extracting %d data frames before and after probe\n", o->padFrames);
    
    int subject = nes_SUBJECT_ID; /* Only one in pre-test. */
    fprintf(stderr, "output files are named \"%s/s%03d_r{BBB}{NN}.txt\"\n", o->outDir, subject);

    int nc = nes_NUM_INPUT_CHANNELS;       /* Total count of data channels (incl. SBS). */
    int ne = nes_NUM_SCALP_ELECTRODES;     /* Count of real electrode channels. */
    assert(ne == nc - 2);                  /* One synthetic channel and the SBS signal. */
    int ic_syn = nc - 2;                   /* Index of synthetic channel. */
    int ic_sbs = nc - 1;                   /* Index of SBS channel. */
    int nt = nes_NUM_INPUT_FRAMES;         /* Count of data frames expected in input file. */
    int np = nes_NUM_PHASES_PER_SESSION;   /* Number of PHASES per file. */
    int nb = nes_NUM_BLOCKS_PER_PHASE;     /* Number of blocks pr phase. */
    int nr = nes_NUM_PROBES_PER_BLOCK;     /* Number of probes per block. */ 
    
    double fsmp = nes_SAMPLING_FREQ;
    
    /* Read the input files: */
    int32_t nt_skip = 0; /* Number of data frames to skip. */
    int32_t nt_take = 0; /* To real all lines. */
    int32_t nl_read = 0; /* Input lines read. */
    int32_t nt_read = 0; /* Input data frames read. */
    double **smp = neuromat_eeg_data_read(stdin, nt_skip, nt_take, nc, &nl_read, &nt_read);
    fprintf(stderr, "read %d lines, %d frames.\n", nl_read, nt_read);
    assert(nt_read == nt);
    
    /* Find the PROBE segments in the file: */
    int kt_ini[np*nb*nr], kt_fin[np*nb*nr]; /* Initial and final frame index of each PROBE. */
    nes_get_PROBES(nt, nc, smp, fsmp, ic_sbs, np, nb, nr, kt_ini, kt_fin);
    
    /* Process and write the data: */
    int32_t nt_PROBE_pad = o->padFrames;             /* Number of padding frames to add per PROBE. */
    int32_t nt_PHASE_pad = 2 * o->padFrames;         /* Number of padding frames to add per PHASE. */
    int32_t nt_PROBE_dat = nes_NUM_FRAMES_PER_PROBE; /* Number of data (stimulus) frames per PROBE. */
    for (int32_t ip = 0; ip < np; ip++)
      { int32_t it_ini = kt_ini[ip*nb*nr] - nt_PHASE_pad; /* Initial frame of PHASE. */
        int32_t it_fin = kt_fin[ip*nb*nr] + nt_PHASE_pad; /* Final frame of PHASE. */
        assert((it_ini >= 0) && (it_fin < nt));
        
        /* Filter phase and write it out for debugging: */
        nes_PHASE_filter(nt, nc, smp, fsmp, ic_sbs, ip, it_ini, it_fin);
        nes_PHASE_write(o->outDir, nt, nc, smp, fsmp, ic_syn, ic_sbs, subject, ip, it_ini, it_fin);
        
        /* Split out PROBE files: */
        for (int32_t ib = 0; ib < nb; ib++)
          { for (int32_t ir = 0; ir < nr; ir++)
              { int32_t cr = ir + (ib + ip*nb)*nr; /* Index of PROBE in SESSION. */
                int32_t ntr = kt_fin[cr] + 1 - kt_ini[cr]; /* Num frames in PROBE according to SBS. */
                assert((ntr >= nt_PROBE_dat - 25) && (ntr <= nt_PROBE_dat + 100)); /* Real frame count cannot be too far from expected. */
                int32_t jt_ini = kt_ini[cr] - nt_PROBE_pad;
                int32_t jt_fin = kt_ini[cr] + nt_PROBE_dat + nt_PROBE_pad; /* Cut probe to nominal size, not SBS size. */
                nes_PROBE_write(o->outDir, nt, nc, smp, fsmp, ic_syn, ic_sbs, subject, ip, ib, ir, jt_ini, jt_fin);
              }
          }
      }
    return 0;
  }
    
void nes_get_PROBES
  ( int32_t nt, 
    int32_t nc, 
    double **smp, 
    double fsmp,
    int32_t ic_sbs,
    int32_t np, 
    int32_t nb, 
    int32_t nr,
    int32_t kt_ini[], /* (OUT) */
    int32_t kt_fin[]  /* (OUT) */
  )
  {
    assert((0 < ic_sbs) && (ic_sbs < nc));
    
    int32_t kt = 0; /* Index of next frame to be parsed. */
    int32_t cr = 0; /* Total count of PROBEs found so far: */
    int32_t nt_3s = (int32_t)floor(3*fsmp + 0.5); /* Three seconds of frames. */
    int32_t nt_30s = (int32_t)floor(30*fsmp + 0.5); /* Thirty seconds of frames. */
    for (int32_t ip = 0; ip < np; ip++)
      { for (int32_t ib = 0; ib < nb; ib++)
          { for (int32_t ir = 0; ir < nr; ir++) 
              { /* Search for start of next PROBE: */
                while ((kt < nt) && (smp[kt][ic_sbs] == 0.5)) { kt++; }
                assert (kt < nt);
                assert ((smp[kt][ic_sbs] == 0.0) || (smp[kt][ic_sbs] == 1.0));
                kt_ini[cr] = kt;
                /* Search for end of PROBE: */
                while ((kt < nt) && (smp[kt][ic_sbs] != 0.5)) { kt++; }
                assert (kt < nt);
                assert (smp[kt][ic_sbs] == 0.5);
                kt_fin[cr] = kt - 1;
                int32_t ntg = kt_ini[cr] - ( cr == 0 ? 0 : kt_fin[cr-1] + 1 ); /* Frames in SETUP/REST/GAP before this PROBE. */
                int32_t ntr = kt_fin[cr] + 1 - kt_ini[cr]; /* Number of frames in probe. */
                fprintf(stderr, "PROBE %4d = p%d_b%d_r%02d : %8d fill %8d frames %8d ..%8d\n", cr, ip, ib, ir, ntg, ntr, kt_ini[cr], kt_fin[cr]);
                /* Consistency checks: */
                if (ir == 0)
                  { /* Previous interval must be SETUP or REST, at least ~30 seconds: */
                    if (cr > 0) { assert(ntr >= nt_30s - 25); }
                  }
                else
                  { /* Previous interval must be at least ~3 seconds: */
                    if (cr > 0) { assert(ntg >= nt_3s - 25); }
                    /* Number of frames must be about 30 seconds: */
                    assert ((ntr >= nt_30s - 25) && (ntr <= nt_30s + 100));
                  }
                /* Found one more run: */
                cr++;
              }
          }
      }
    int32_t ntg = (nt + 1) - ( cr == 0 ? 0 : kt_fin[cr-1] + 1 ); /* Frames in SETUP/REST/GAP after last PROBE. */
    fprintf(stderr, "%8d fill at end of file\n", ntg);
    assert(cr == np*nb*nr);
    return;
  }
                
void nes_PHASE_filter
  ( int32_t nt, 
    int32_t nc, 
    double **smp, 
    double fsmp,
    int32_t ic_sbs,
    int32_t ip,
    int32_t it_ini,
    int32_t it_fin
  )
  {
    assert((0 < ic_sbs) && (ic_sbs < nc));
    assert((it_ini >= 0) && (it_fin < nt));
    
    int32_t ns = it_fin + 1 - it_ini; /* Number of frames in PHASE. */
    
    /* Allocate the FFT work areas and plans: */
    double *x_time = (double*) fftw_malloc(sizeof(double)*(ns+4)); /* Time-domain signal. */
    double *x_freq = (double*) fftw_malloc(sizeof(double)*(ns+4)); /* Hartley transform. */
    fftw_plan pdir = fftw_plan_r2r_1d(ns, x_time, x_freq, FFTW_DHT, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    fftw_plan pinv = fftw_plan_r2r_1d(ns, x_freq, x_time, FFTW_DHT, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

    double *F = nes_make_filter(ns, fsmp); /* Fourier domain filter. */
    
    for (int32_t ic = 0; ic < nc; ic++) 
      { if (ic != ic_sbs)
          { /* Extract the data to be filtered: */
            for (int32_t it = 0; it < ns; it++) { x_time[it] = smp[it_ini + it][ic]; }

            /* Apply filter: */
            nes_signal_filter(ns, fsmp, x_time, x_freq, F, pdir, pinv);

            /* Return filtered signal & trend to {smp}, scaled to preserve total power: */
            for (int32_t it = 0; it < ns; it++) { smp[it_ini + it][ic] = x_time[it]/nt; }
          }
      }
      
    fftw_destroy_plan(pdir);
    fftw_destroy_plan(pinv);
    fftw_free(x_freq);
    fftw_free(x_time);
    free(F);
    return;
  }
    
void nes_signal_filter(int32_t ns, double fsmp, double x[], double X[], double F[], fftw_plan pdir, fftw_plan pinv)
  {
    /* Convert to frequency domain: */
    fftw_execute(pdir);

    /* Apply frequency filter: */
    int kf0;
    for (kf0 = 0; kf0 <= ns-kf0; kf0++) 
      { int kf1 = (ns - kf0) % ns; /* The other Hartley element with same absolute freq. */
        assert(F[kf1] == F[kf0]); /* Filter must be even. */
        double Fk = F[kf0];
        X[kf0] = X[kf0]*Fk;
        if (kf1 != kf0) { X[kf1] = X[kf1]*Fk; }
      }

    /* Return to time domain: */
    fftw_execute(pinv);
  }

void nes_PHASE_write
  ( char *outDir,
    int32_t nt, 
    int32_t nc, 
    double **smp, 
    double fsmp,
    int32_t ic_syn,
    int32_t ic_sbs,
    int32_t subject,
    int32_t ip,
    int32_t it_ini,
    int32_t it_fin
  )
  {
    assert((0 < ic_sbs) && (ic_sbs < nc));
    assert((it_ini >= 0) && (it_fin < nt));
    
    int32_t ntp = it_fin + 1 - it_ini; /* Number of frames in PHASE. */

    char *fname = NULL;
    asprintf(&fname, "%s/s%03d_p%d.txt", outDir, subject, ip);
    FILE *wr = open_write(fname, TRUE);
    
    int32_t it_step = 1; /* Subsampling step. */

    char *ptype = ( ip == 0 ? "DBS0" : "DBS1" );
    nes_header_write(wr, nt, nc, fsmp, it_step, ic_syn, ic_sbs, subject, ptype, it_ini, it_fin);

    for (int32_t it = 0; it < ntp; it = it + it_step)
      { double *vi = smp[it_ini + it];
        for (int32_t ic = 0; ic < nc; ic++)
          { fprintf(wr, "%9.3f", vi[ic]);
            if (ic > 0) { fputc(' ', wr); }
          }
        fputc('\n', wr);
      }
    fclose(wr);
    free(fname);
  }
        
void nes_PROBE_write
  ( char *outDir,
    int32_t nt, 
    int32_t nc, 
    double **smp, 
    double fsmp,
    int32_t ic_syn,
    int32_t ic_sbs,
    int32_t subject,
    int32_t ip,
    int32_t ib,
    int32_t ir,
    int32_t it_ini,
    int32_t it_fin
  )
  {
    assert((0 < ic_sbs) && (ic_sbs < nc));
    assert((it_ini >= 0) && (it_fin < nt));

    int32_t ntr = it_fin + 1 - it_ini; /* Number of frames in PHASE. */
    
    int32_t itm = (it_ini + it_fin)/2; /* A frame in the middle of the run. */
    assert((smp[itm][ic_sbs] == 0.0) || (smp[itm][ic_sbs] == 1.0));

    char *fname = NULL;
    asprintf(&fname, "%s/s%03d_p%d_b%d_r%02d.txt", outDir, subject, ip, ib+1, ir+1);
    FILE *wr = open_write(fname, TRUE);
    
    int32_t it_step = 4; /* Subsampling step. */

    char *ptype = ( ip == 0 ? "DBS0" : "DBS1");
    char *rtype = ( smp[itm][ic_sbs] == 0 ? "HS" : "HD" );
    char *type = txtcat(ptype, rtype);
    nes_header_write(wr, nt, nc, fsmp, it_step, -1, ic_sbs, subject, type, it_ini, it_fin);

    for (int32_t it = 0; it < ntr; it = it + it_step)
      { double *vi = smp[it_ini + it];
        for (int32_t ic = 0; ic < nc; ic++)
          { if (ic != ic_syn)
              { fprintf(wr, "%9.3f", vi[ic]);
                if (ic > 0) { fputc(' ', wr); }
              }
          }
        fputc('\n', wr);
      }
    fclose(wr);
    free(fname);
  }
        
void nes_header_write
  ( FILE *wr, 
    int32_t nt, 
    int32_t nc, 
    double fsmp,
    int32_t it_step,
    int32_t ic_syn,
    int32_t ic_sbs,
    int32_t subject,
    char *type,
    int32_t it_ini,
    int32_t it_fin
  )
  {
    assert((0 < ic_sbs) && (ic_sbs < nc));
    assert((it_ini >= 0) && (it_fin < nt));
    
    int32_t ns_orig = it_fin + 1 - it_ini; /* Number of frames in range, before downsampling. */
    int32_t ns_file = ns_orig/it_step; /* Number of samples really written. */

    assert(nc == 21);
    int32_t ne = 19;
    
    int ncw = (ic_syn >= 0 ? nc : nc-1); /* Number of channels in output file. */
    
    neuromat_eeg_header_t *hot = neuromat_eeg_header_new();
    hot->nt = ns_file;
    hot->nc = 0;  /* To be updated below. */
    hot->ne = 0;  /* To be updated below. */
    hot->fsmp = fsmp/it_step;
    hot->type = type;
    hot->subject = subject;
    
    hot->orig = neuromat_eeg_source_new();
    hot->orig->nt = nt;
    hot->orig->it_ini = it_ini;
    hot->orig->it_fin = it_fin;
    hot->orig->fsmp = fsmp;
    hot->orig->subject = subject;
    hot->orig->run = INT32_MIN; /* All runs. */
    
    assert(ne == 19);
    char *ename[] = { "Fp1", "Fp2", "F3", "F4", "C3", "C4", "P3", "P4", "O1", "O2", "F7", "F8", "T3", "T4", "T5", "T6", "Fz", "Cz", "Pz" }; 
    int32_t res;
    for (int32_t ie = 0; ie < ne; ie++) { res = neuromat_eeg_header_append_electrode_channel(hot, ename[ie]); assert (res == ie); }
    if (ic_syn >= 0) { res = neuromat_eeg_header_append_electrode_channel(hot, "SYN"); assert(res == ic_syn); }
    res = neuromat_eeg_header_append_marker_channel(hot, "SBS"); assert(res == ncw - 1);

    fprintf(stderr, "writing header...\n");
    neuromat_eeg_header_write(wr, hot);
    /* neuromat_eeg_header_free(hot); */
  }

double *nes_make_filter(int32_t nt, double fsmp)
  {
    /* Bandpass filter parameters: */
    double f0_kill = 1.0;  /* Low kill frequency. */
    double f0_pass = 2.0;  /* Low pass frequency. */
    double f1_pass = 45.0; /* High pass frequency (Hz). */
    double f1_kill = 60.0; /* High kill frequency (Hz). */
    
    /* Notch filter components: */
    double nn = 3; /* Number of notch filters to use. */
    double fn[] = { 20.1465, 29.3009, 40.2935 };  /* Center frequencies (Hz). */
    double fd[] = {  0.0900,  0.0200,  0.1800 };  /* Relative widths. */
    
    double *F = notnull(malloc(nt*sizeof(double)), "no mem");
    int32_t kf_max = nt/2; /* Max digital frequency in file. */
    for (int32_t kf0 = 0; kf0 <= kf_max; kf0++)
      { double f = kf0*fsmp/nt; /* Frequency of Hartley component in Hz. */
        double F_hi = 1.0 - neuromat_filter_lowpass_sigmoid(f, f0_kill, f0_pass); /* Highpass. */
        double F_lo = neuromat_filter_lowpass_sigmoid(f, f1_pass, f1_kill);     /* Lowpass. */
        double Fk = F_lo * F_hi;
        /* Apply the notch filters: */
        for (int32_t in = 0; in <nn; in++)
          { double dz = fabs((f - fn[in])/fd[in]);
            double F_ni = (dz > 8.0 ? 0.0 : 1.0 - exp(-0.5*dz*dz));
            Fk = Fk * F_ni;
          }
        F[kf0] = Fk;
        int32_t kf1 = (nt - kf0) % nt; /* Conjugated digital frequency. */
        if (kf1 != kf0) { F[kf1] = Fk; }
      }
    return F;
  }
#define nes_MAX_BLOCK_ID (nes_NUM_BLOCKS_PER_PHASE)
  /* Max BLOCK index. */

#define nes_MAX_PROBE_ID (nes_NUM_PROBES_PER_BLOCK)
  /* Max PROBE index. */

#define nes_MAX_PAD_FRAMES (1250)
  /* Max allowed frame pading (2.5 seconds, to fit in 3 seconds of GAP). */

nes_options_t *nes_parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    nes_options_t *o = notnull(malloc(sizeof(nes_options_t)), "no mem"); 
    
    /* Parse keyword parameters: */
    
    if (argparser_keyword_present(pp, "-firstBlockID"))
      { o->firstBlockID = (int)argparser_get_next_int(pp, 1, nes_MAX_BLOCK_ID); }
    else
      { o->firstBlockID = 1; }
      
    if (argparser_keyword_present(pp, "-firstProbeID"))
      { o->firstProbeID = (int)argparser_get_next_int(pp, 1, nes_MAX_PROBE_ID); }
    else
      { o->firstProbeID = 1; }
      
    /* Parse trigger channel name option and look it up in channel name list: */
    argparser_get_keyword(pp, "-padFrames");
    o->padFrames = (int)argparser_get_next_int(pp, 1, nes_MAX_PAD_FRAMES);
    
    argparser_get_keyword(pp, "-outDir");
    o->outDir = argparser_get_next_non_keyword(pp);
  
    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }
