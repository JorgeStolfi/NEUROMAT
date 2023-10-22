#define PROG_NAME "nmeeg_filter"
#define PROG_DESC "Applies a bandpass filter to an EEG dataset."
#define PROG_VERS "2013-11-15"

#define nmeeg_filter_C_COPYRIGHT \
  "Copyright © 2013 by the State University of Campinas (UNICAMP)"

/* Last edited on 2023-10-21 21:56:18 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    [ -filter {F_TYPE} {F_LO0} {F_LO1} {F_HI1} {F_HI0} ] \\\n" \
  "    [ -invert {INV_FLAG} ] \\\n" \
  "    [ -resample {R_STEP} ] \\\n" \
  "    [ -trend {TREND_DEG} {TREND_KEEP} ] \\\n" \
  "    [ -verbose ] \\\n" \
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
  "  The program reads from standard input an EEG dataset in the plaintext (\".txt\") format," \
  " applies a polynomial-fit and/or frequency bandpass/bandkill filter to each electrode channel, and writes the" \
  " result (possibly sub-sampled) to standard output.\n" \
  "\n" \
  "  The input file should contain a certain" \
  " number {NT} of data frames, each with the same number {NC} of channel" \
  " values, including {NE} electrode readings and {NC-NE} non-electrode" \
  " channels (such as trigger and phase markers).\n" \
  "\n" \
  "  The input file must have header records that specify, among other things, the number of" \
  " channels {NC} and their names {CHNAME[0]} {CHNAME[1]} ... {CHNAME[NC-1]}," \
  " the number of data frames {NT}, the sampling frequency {FSMP}, and" \
  " the number of electrode channels {NE}.  The first {NE} channels" \
  " in each data frame are assumed to be electrode" \
  " potentials.\n" \
  "\n" \
  "  The output file will have the same format, except that the number" \
  " of frames may be reduced and the header will include information about the filter band.\n" \
  "\n" \
  "TREND PROCESSING\n" \
  "  The frequency bandpass or bandkill filters used in this program assume that the" \
  " signal is periodic, and that the input samples span one full period.  When" \
  " such filters are used on a finite segment of a non-periodic signal, artifacts" \
  " usually appear at each end due to the implicit jump from the last to the" \
  " first sample.\n" \
  "\n" \
  "  To reduced such artifacts, the program can be asked to fit a" \
  " `trend' polynomial {P} to each channel, subtract {P} from the channel" \
  " before filtering, and optionally add {P} back after filtering (but before" \
  " downsampling, if any).  See the \"-trend\" option below.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -filter {F_TYPE} {F_LO0} {F_LO1} {F_HI1} {F_HI0}\n" \
  "    This optional argument specifies the frequency bandpass/bandkill filter type {F_TYPE} and" \
  " the band parameters {F_LO0} {F_LO1} {F_HI1} {F_HI0}.\n" \
  "\n" \
  "    The filter's gain is always real and non-negative for all frequencies.  Roughly, in" \
  " bandpass mode, the filter's gain will be very small for frequencies" \
  " below {F_LO0} or above {F_HI0}, and nearly 1 for frequencies" \
  " between {F_LO1} and {F_HI1}.  The transitions between {F_LO0} and {F_LO1} and" \
  " between {F_HI0} and {F_HI1} depend on the filter type.  In bandkill mode, the gain is complemented" \
  " relative to unit gain.\n" \
  "\n" \
  "    All frequencies are" \
  " in hertz (cycles per second).   If this option is not specified, the" \
  " gain is 1 for all frequencies in bandpass mode, and 0 for all" \
  " frequencies in bandkill mode.\n" \
  "\n" \
  "    Currently defined filter types are \"G\" (difference" \
  " of Gaussians, recommended), \"SG\" (sigmoid transitions" \
  " in log frequency domain), \"BUR8\" (real Butterworth filter" \
  " of order 8), and \"BUC12\" (complex Butterworth filter" \
  " of order 12).  The type \"D\" is also valid and specifies" \
  " a delay filter, useful for testing.\n" \
  "\n" \
  "  -invert {INV_FLAG}\n" \
  "    This optional argument specifies whether the filter is to" \
  " be used in the natural `bandpass' mode (if {INV_FLAG} is 0 or \"F\") or" \
  " in `bandkill' mode (if {INV_FLAG} is 1 or \"T\").  In the second" \
  " case, the gain for each frequency is defined as 1 minus the" \
  " natural bandpass gain.  If omitted, the program assumes bandpass ({INV_FLAG=0}).\n" \
  "\n" \
  "  -resample {R_STEP}\n" \
  "    This optional argument specifies that the output should be" \
  " downsampled after filtering, by taking one every {R_STEP} samples.  If" \
  " omitted, {R_STEP=1} is assumed (no downsampling).\n" \
  "\n" \
  "  -trend {TREND_DEG} {TREND_KEEP}\n" \
  "    This optional argument specifies whether a polynomial trend should be" \
  " separately handled for each electrode signal.  If {TREND_DEG} is" \
  " zero or positive, the program will fit a polynomial of of" \
  " degree {TREND_DEG} to each electrode signal, and subtract it from from the" \
  " signal, before applying the Fourier filter.  In that case, the second" \
  " argument {TREND_KEEP} specifies whether the fitted polynomial trend" \
  " is to be added back to the filtered signal (\"T\" or 1), or" \
  " discarded (\"F\" or 0).  If {TREND_DEG} is [-1} (or if the option is not" \
  " present), this trend fittting and separation process is not performed, and" \
  " the {TREND_KEEP} argument is ignored.\n" \
  "\n" \
  "  -verbose\n" \
  "    This optional flag causes some additional diagnostic output to be printed.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  neuromat_eeg_plot_signals.sh(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2013-06-03 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2013-11-15 Converted to {argparser.h}, added \"-pad\" flag.\n" \
  "  2013-11-18 Replaced the \"-pad\" flag by \"-trend\" flag.\n" \
  "  2013-11-20 Added keep/remove option to \"-trend\" flag.\n" \
  "  2013-11-20 Made \"-filter\" optional.\n" \
  "  2013-11-21 Added \"-verbose\".\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " nmeeg_filter_C_COPYRIGHT ".\n" \
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

#include <r2.h>
#include <affirm.h>
#include <jsfile.h>
#include <argparser.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_io.h>
#include <neuromat_eeg_header.h>
#include <neuromat_filter.h>
#include <neuromat_eeg_stats.h>

typedef struct nef_freq_filter_t
  { char *type;
    double flo0;
    double flo1;
    double fhi1;
    double fhi0;
  } nef_freq_filter_t;
  /* Type and characteristic frequencies of the frequency filter. If a filter was 
    not specified, the band paramters are {0,0,+INF,+INF}, and {type} is NULL. */

typedef struct nef_options_t
  { int trend_deg;             /* Degree of polynomial trend to separate, or {-1} for no trend. */
    bool_t trend_keep;         /* Tells whether trend polynomial is to be kept (true) of discarded (false). */
    nef_freq_filter_t filter;  /* Freq filter parameters. */
    bool_t invert;             /* Freq filter mode: FALSE for bandpass, TRUE for bandkill. */
    int resample;              /* Resampling step. */
    bool_t verbose;            /* True to print more info. */
  } nef_options_t;
  /* Arguments from command line. */

void nef_write_eeg_signals
  ( int nt, 
    int nc,
    double **val, 
    int ne,
    double fsmp,
    int it_ini,
    int it_fin,
    int it_step,
    neuromat_eeg_header_t *h
  );
  /* Writes the signals {val[0..nt-1][0..nc-1]}, comprising {nc}
    channels sampled at {nt} times, to standard output.
    Writes only one very {it_step} frames, beginning with frame
    {val[it_ini]} and ending with frame {val[it_fin]}.
    
    Writes the header {h} to the file, setting its fields
    {h->nt,h->nc,h->ne} to {nw,nc,ne} where {nw} is the number of
    subsampled frames. Also set {h->kfmax} to a null value. Assumes that
    {fsmp} refers to the full dataset {val[0..nt-1]}, so sets {h->fsmp}
    to {fsmp/it_step}. Assumes that {h->orig->it_ini} and {h->orig->it_fin}
    describe the position of {val[0..nt-1]} in the original file, so
    adjusts those fields to account for the clipping and subsampling
    defined by {it_ini,it_fin}. */
  
void nef_update_trend_filter_in_header(neuromat_eeg_header_t *h, int tdeg, bool_t tkeep);
  /* Updates the trend filter parameters in the header {h} to account for the
    trend processing defined by the parameters {tdeg, tkeep}. */

void nef_update_freq_filter_in_header
  ( neuromat_eeg_header_t *h, 
    double flo0, 
    double flo1, 
    double fhi1,
    double fhi0,
    int finvert
  );
  /* Updates the freq filter parameters in the header {h} to account for the
    soft bandpass/bandkill filtering with parameters {flo0,flo1,fhi1,fhi0,finvert}.
    Expects the input data to be unfiltered. If the header already has
    filter parameters set, tries to merge them, but the result may not
    correctly describe the effect of the two sucessive filtering steps. */
      
nef_options_t *nef_parse_options(int argc, char **argv);
  /* Parses the command line arguments. */
  
int main(int argc, char **argv);
  /* Main prog. */

#define nef_MAX_RESAMPLE 10
  /* Max downsampling step allowed (paranoia). */

#define nef_MAX_TREND_DEG 5
  /* Max downsampling step allowed (paranoia). */

int main(int argc, char **argv)
  {
    nef_options_t *o = nef_parse_options(argc, argv);
    
    nef_freq_filter_t *ff = &(o->filter);
    
    if (ff->type == NULL)
      { fprintf(stderr, "no frequency filtering\n"); }
    else
      { fprintf(stderr, "freq filter type %s", ff->type);
        fprintf(stderr, " low freq ramp [%.6f _ %.6f]", ff->flo0, ff->flo1);
        fprintf(stderr, " high freq ramp [%.6f _ %.6f]\n", ff->fhi1, ff->fhi0);
        fprintf(stderr, "filter mode is %s\n", (o->invert ? "bandpass" : "bandkill"));
      }

    /* Resampling: */
    fprintf(stderr, "output will have one every %d samples\n", o->resample);
    demand(o->resample > 0, "invalid resampling step");
    
    fprintf(stderr, "writing result to standard output\n");

    /* Read the EEG header, provide defaults: */
    int nl = 0;
    neuromat_eeg_header_t *h = neuromat_eeg_header_read(stdin, 20, 600.0, &nl);
    fprintf(stderr, "read %d header lines\n", nl);
    
    /* Get important parameters: */
    int nc = h->nc; /* Number of channels per data frame (including triggers etc.). */
    int ne = h->ne; /* Assumes that the first {ne} channels are electrode potentials. */
    double fsmp = h->fsmp;
    char **chname = h->chname;
    fprintf(stderr, "input file has %d channels, including %d electrode potentials\n", nc, ne);
    fprintf(stderr, "input sampling frequency %.10g Hz\n", fsmp);

    /* Read the EEG data frames: */
    int nl0 = nl;
    int nt = 0;
    double **val = neuromat_eeg_data_read(stdin, 0, 0, nc, &nl, &nt);
    fprintf(stderr, "read %d lines, got %d data frames\n", nl - nl0, nt);
    demand(nt > 0, "no frames read");
    demand(nt == h->nt, "inconsistent header {nt}");
    
    /* Filter procedures: */        
    auto complex filter_gain_unit(int kf, int nf, double fsmp);
    auto complex filter_gain_butterworth8(int kf, int nf, double fsmp);
    auto complex filter_gain_cbutterworth12(int kf, int nf, double fsmp);
    auto complex filter_gain_delay(int kf, int nf, double fsmp);
    auto complex filter_gain_gaussian(int kf, int nf, double fsmp);
    auto complex filter_gain_biquadratic(int kf, int nf, double fsmp);
    auto complex filter_gain_bisigmoid(int kf, int nf, double fsmp);
  
    /* Choose filter procedure and max significant freq {fcut} on bandpass output: */
    neuromat_filter_t *filter_gain = NULL;
    double fcut = NAN;
    
    if (ff->type == NULL)
      { /* No frequency filter: */
        filter_gain = &filter_gain_unit;
        fcut = +INF;
      }
    else if (strcmp(ff->type, "SG") == 0)
      { /* A bisigmoid {flo0,flo1,fhi1,fhi0}: */
        filter_gain = &filter_gain_bisigmoid;
        fcut = ff->fhi0;
      }
    else if (strcmp(ff->type, "DL") == 0)
      { /* The 1 sec delay filter for testing: */
        filter_gain = &filter_gain_delay;
        fcut = +INF;
      }
    else if (strcmp(ff->type, "G") == 0)
      { /* The double Gaussian filter defined by {flo0,fhi0}: */
        filter_gain = &filter_gain_gaussian;
        fcut = ff->fhi0;
      }
    else if (strcmp(ff->type, "BQ") == 0)
      { /* The double biquadratic {fsmp} filter: */
        filter_gain = &filter_gain_biquadratic;
        fcut = +INF;
      }
    else if (strcmp(ff->type, "BUR8") == 0)
      { /* The real Butterworth filter of order 8 with cutoff {~1.2*fhi1}: */
        filter_gain = &filter_gain_butterworth8;
        fcut = 2*ff->fhi1;
      }
    else if (strcmp(ff->type, "BUC12") == 0)
      { /* The complex Butterworth filter of order 12: */
        filter_gain = &filter_gain_cbutterworth12;
        fcut = 2*ff->fhi1;
      }
    else
      { demand(FALSE, "unknown filter type"); }
    /* In any case the max significant freq is the Nyquist limit: */
    fcut = fmin(fcut, fsmp/2);
      
    /* Compute subsampling indexing parameters: */
    if ((o->invert != 0) || (fcut >= fsmp/2))
      { demand(o->resample == 1, "no resampling allowed with bandkill filter"); }
    else 
      { int max_resample = (int)floor((fsmp/2)/fcut); 
        if (max_resample > nef_MAX_RESAMPLE) { max_resample = nef_MAX_RESAMPLE; }
        demand(o->resample < max_resample, "resampling step too large for this filter");
      }
    
    /* Apply the filter: */
    neuromat_filter_apply(nt, ne, val, fsmp, o->trend_deg, o->trend_keep, filter_gain, o->verbose);

    /* Remember the filtering parameters in header: */
    nef_update_trend_filter_in_header(h, o->trend_deg, o->trend_keep);
    nef_update_freq_filter_in_header(h, ff->flo0, ff->flo1, ff->fhi1, ff->fhi0, o->invert);

    /* Subsampling parameters: */
    int it_step = o->resample;
    int it_ini = it_step/2;
    int nw = (nt + it_step - 1 - it_ini)/it_step; /* Number of frames to write. */
    int it_fin = it_ini + it_step*((nt - it_ini - 1)/it_step);
    assert((it_fin >= nt - it_step) && (it_fin < nt));
    fprintf(stderr, "saving %d frames (%d to %d step %d)\n", nw, it_ini, it_fin, it_step);

    /* Compute basic statistics of filtered electrode signals: */
    double vavg[nc], vvar[nc], vmin[nc], vmax[nc];
    double vdev_max, vmin_min, vmax_max;  /* Maximum variance among electrode channels. */
    bool_t zeroMean = FALSE;
    neuromat_eeg_stats_per_channel(nt, nc, val, zeroMean, vavg, vvar, vmin, vmax, ne, &vmin_min, &vmax_max, &vdev_max);
    fprintf(stderr, "--- filtered signal statistics ---\n");
    neuromat_eeg_stats_per_channel_print(stderr, nc, chname, vavg, vvar, vmin, vmax);
    fprintf(stderr, "range of electrode values = [ %8.5f _ %8.5f ]\n", vmin_min, vmax_max);
    fprintf(stderr, "maximum dev = %8.5f\n", vdev_max); 
    fprintf(stderr, "\n");

    /* Write them out, subsampled: */
    nef_write_eeg_signals(nt, nc, val, ne, fsmp, it_ini, it_fin, it_step, h);
      
    int it;
    for (it = 0; it < nt; it++) { free(val[it]); } 
    free(val);
    
    return 0;
    
    complex filter_gain_unit(int kf, int nf, double fsmp)
      { 
        double w = 1.0;
        assert(ff->flo0 == 0);
        assert(ff->flo1 == 0);
        assert(ff->fhi1 == +INF);
        assert(ff->fhi0 == +INF);
        return (o->invert == 0 ? w : 1-w);
      }
    
    complex filter_gain_bisigmoid(int kf, int nf, double fsmp)
      { 
        demand((kf >= 0) && (kf < nf), "bad {kf}");
        if (2*kf > nf) { kf = nf - kf; } /* Reduce {kf} to absolute freq in {0..nf/2} */
        double f = kf*fsmp/nf;
        double wlo = neuromat_filter_lowpass_sigmoid(f, ff->flo0, ff->flo1);
        double whi = neuromat_filter_lowpass_sigmoid(f, ff->fhi1, ff->fhi0);
        double w = (1-wlo)*whi;
        return (o->invert == 0 ? w : 1-w);
      }
    
    complex filter_gain_butterworth8(int kf, int nf, double fsmp)
      { 
        double btfhi = ff->fhi1/0.8; btfhi = floor(btfhi + 0.5);   /* Highpass shoulder. */
        double btflo = ff->flo1*0.8; btflo = floor(btflo*10.0 + 0.5)/10.0; /* Lowpass shoulder. */
        int btorder = 8;
        demand((kf >= 0) && (kf < nf), "bad {kf}");
        if (2*kf > nf) { kf = nf - kf; } /* Reduce {kf} to absolute freq in {0..nf/2} */
        double f = kf*fsmp/nf;
        double wlo = neuromat_filter_lowpass_butterworth(f, btflo, btorder);
        double whi = neuromat_filter_lowpass_butterworth(f, btfhi, btorder);
        double w = (1-wlo)*whi;
        return (o->invert == 0 ? w : 1-w);
      }
    
    complex filter_gain_cbutterworth12(int kf, int nf, double fsmp)
      { 
        double btfhi = ff->fhi1/0.8; btfhi = floor(btfhi + 0.5);   /* Highpass shoulder. */
        double btflo = ff->flo1*0.8; btflo = floor(btflo*10.0 + 0.5)/10.0; /* Lowpass shoulder. */
        int btorder = 12;
        demand((kf >= 0) && (kf < nf), "bad {kf}");
        if (2*kf > nf) { kf = kf - nf; } /* Reduce {kf} to {-nf/2..nf/2} */
        double f = kf*fsmp/nf;
        complex wlo = neuromat_filter_lowpass_cbutterworth(f, btflo, btorder);
        complex whi = neuromat_filter_lowpass_cbutterworth(f, btfhi, btorder);
        complex w = (1-wlo)*whi;
        /* if (kf == 0) { fprintf(stderr, "  %10s  %25s %25s %25s\n", "freq", "wlo", "whi", "w"); } */
        /* fprintf(stderr, "  %10.3f", f); */
        /* fprintf(stderr, "  %12.8f %12.8f", creal(wlo), cimag(wlo)); */
        /* fprintf(stderr, "  %12.8f %12.8f", creal(whi), cimag(whi)); */
        /* fprintf(stderr, "  %12.8f %12.8f", creal(w), cimag(w)); */
        /* fprintf(stderr, "\n"); */
        return (o->invert == 0 ? w : 1-w);
      }
    
    complex filter_gain_delay(int kf, int nf, double fsmp)
      { 
        demand((kf >= 0) && (kf < nf), "bad {kf}");
        complex w;
        if (kf == 0)
          { w = 1; }
        else if (2*kf == nf)
          { w = 0; }
        else 
          { if (2*kf > nf) { kf = kf - nf; }
            int kdelay = (int)floor(fsmp);
            if (2*kdelay >= nf) { kdelay = nf/4; }
            w = cexp(-kdelay*kf*2*M_PI/nf*I);
          }
        return (o->invert == 0 ? w : 1-w);
      }

    complex filter_gain_gaussian(int kf, int nf, double fsmp)
      { 
        demand((kf >= 0) && (kf < nf), "bad {kf}");
        double w;
        if (kf == 0) 
          { w = 1; }
        else
          { double f = kf*fsmp/nf;
            double whi = (ff->fhi0 >= +INF ? 1.0 : neuromat_filter_lowpass_gauss(f, ff->fhi0, 0.001, fsmp));
            double wlo = (ff->flo1 <= 0 ? 1 : 1 - neuromat_filter_lowpass_gauss(f, ff->flo1, 0.001, fsmp));
            w = whi*wlo;
          }
        return (o->invert == 0 ? w : 1-w);
      }        

    complex filter_gain_biquadratic(int kf, int nf, double fsmp)
      { 
        demand((kf >= 0) && (kf < nf), "bad {kf}");
        double w;
        if (kf == 0) 
          { w = 1; }
        else
          { double f = kf*fsmp/nf;
            w = neuromat_filter_lowpass_biquadratic(f, fsmp/2);
          }
        return (o->invert == 0 ? w : 1-w);
      }        
  }

void nef_write_eeg_signals
  ( int nt, 
    int nc, 
    double **val, 
    int ne,
    double fsmp,
    int it_ini,
    int it_fin,
    int it_step,
    neuromat_eeg_header_t *h
  )
  {
    demand((it_fin - it_ini) % it_step == 0, "invalid subsampling range/step");
    int nw = (it_fin - it_ini)/it_step + 1;
    
    /* Set the dataset size in header: */
    h->nt = nw;
    h->nc = nc;
    h->ne = ne;
    
    /* Adjust original dataframe index range in header to account for subsampling: */
    h->orig->it_ini += it_ini;
    h->orig->it_fin = h->orig->it_ini + (it_fin - it_ini);
    
    /* Find the max freq {ftop} present in the input: */
    double ftop = fsmp/2;
    if (! isnan(h->fhi0)) { ftop = fmin(h->fhi0, ftop); }
    
    /* Adjust the sampling frequency to account for subsampling: */
    double fsub = fsmp/it_step;  /* Sampling freq of resampled filtered signal. */
    demand(fsub >= 2*ftop, "resampling freq too low for filter");
    h->fsmp = fsub;
    h->kfmax = INT_MIN; /* Not meaningful for a time-domain dataset. */

    /* Write to {stdout}: */
    FILE *wr = stdout;
    neuromat_eeg_header_write(wr, h);
    neuromat_eeg_data_write(wr, nt, nc, val, it_ini, it_fin, it_step);
    fflush(wr);
  }

void nef_update_trend_filter_in_header(neuromat_eeg_header_t *h, int tdeg, bool_t tkeep)
  {
    if (h->tdeg < 0)
      { /* Input was not trend-filtered: */
        demand(h->tkeep == INT_MIN, "inconsistent trend filter params in header");
        if (tdeg >= 0)
          { /* New trend filter is not trivial, just save it: */
            h->tdeg = tdeg;
            h->tkeep = (tkeep ? 1 : 0);
          }
      }
    else 
      { demand((h->tkeep == 0) || (h->tkeep == 1), "inconsistent trend filter params in header");
        if (tdeg < 0)
          { /* No trend filter applied, keep original: */ }
        else
          { /* Non-trivial filter applied: */
            if ((h->tdeg != tdeg) || (h->tkeep != tkeep))
              { /* !!! Should try to merge !!! */
                demand(FALSE, "incompatible trend filter params");
              }
          }
      }
  }

void nef_update_freq_filter_in_header
  ( neuromat_eeg_header_t *h, 
    double flo0, 
    double flo1, 
    double fhi1, 
    double fhi0, 
    int finvert
  )
  {
    if (h->finvert < 0)
      { /* Input was not freq-filtered. */
        demand(isnan(h->flo0) && isnan(h->flo1) &&isnan(h->fhi1) && isnan(h->fhi0), "inconsistent freq filter params in header");
        if ((finvert == 1) || (flo1 > 0) || (fhi1 < +INF))
          { /* New freq filter is not trivial, just save it: */
            h->finvert = finvert;
            h->flo0 = flo0; h->flo1 = flo1;
            h->fhi1 = fhi1; h->fhi0 = fhi0;
          }
      }
    else 
      { /* Input was freq-filtered. */
        demand((! isnan(h->flo0)) && (! isnan(h->flo1)), "inconsistent low-shoulder freq filter params in header");
        demand((! isnan(h->fhi1)) && (! isnan(h->fhi0)), "inconsistent high-shoulder freq filter params in header");
        demand((h->finvert == 0) || (h->finvert == 1), "invalid freq filter invert flag in header");
        if ((finvert == 0) && (flo1 <= 0) && (fhi1 >= +INF))
          { /* Unit filter applied, keep old filter. */
          }
        else if ((finvert == 1) && (flo1 <= 0) && (fhi1 >= +INF))
          { /* Kill-all filter applied, overwrites old filter: */
            h->finvert = finvert;
            h->flo0 = flo0; 
            h->flo1 = flo1;
            h->fhi1 = fhi1;
            h->fhi0 = fhi0;
          }
        else
          { /* Non-trivial freq filter applied: */
            if (h->finvert != finvert)
              { /* !!! Should try to merge !!! */
                demand(FALSE, "incompatible freq filter inversion flags");
              }
            /* Check for overlapping low-low or high-high shoulders: */
            if ((fmax(h->flo0, flo0) < fmin(h->flo1, flo1))) 
              { fprintf(stderr, "!! warning: overlapping low-freq shoulders"); }
            if ((fmax(h->fhi1, fhi1) < fmin(h->fhi0, fhi0))) 
              { fprintf(stderr, "!! warning: overlapping high-freq shoulders"); }
            /* Intersect the two filters: */
            h->flo0 = fmax(h->flo0, flo0); 
            h->flo1 = fmax(h->flo1, flo1);
            h->fhi1 = fmin(h->fhi1, fhi1);
            h->fhi0 = fmin(h->fhi0, fhi0);
            demand(h->flo0 < h->fhi0, "filtered everything out?");
            demand((h->flo0 <= h->flo1) && (h->fhi1 <= h->fhi0), "something wrong with freq filter params");
            demand((h->flo1 <= h->fhi0) && (h->flo0 <= h->fhi1), "something else wrong with freq filter params");
            /* Check for overlapping low and high shoulders (no flat band): */
            if (h->flo1 > h->fhi1)
              { fprintf(stderr, "!! warning: overlapping high-freq and low-freq shoulders"); }
          }
      }
  }

nef_options_t *nef_parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    nef_options_t *o = notnull(malloc(sizeof(nef_options_t)), "no mem"); 
    
    /* Parse keyword parameters: */
    
    if (argparser_keyword_present(pp, "-filter"))
      { 
        o->filter.type = argparser_get_next(pp);
        o->filter.flo0 = argparser_get_next_double(pp, 0.0, +INF);
        o->filter.flo1 = argparser_get_next_double(pp, o->filter.flo0, +INF);
        o->filter.fhi1 = argparser_get_next_double(pp, o->filter.flo1, +INF);
        o->filter.fhi0 = argparser_get_next_double(pp, o->filter.fhi1, +INF);
      }
    else
      { o->filter.type = NULL;
        o->filter.flo0 = 0.0;
        o->filter.flo1 = 0.0;
        o->filter.fhi1 = +INF;
        o->filter.fhi0 = +INF;
      }
 
    if (argparser_keyword_present(pp, "-invert"))
      { o->invert = argparser_get_next_bool(pp); }
    else
      { o->invert = FALSE; }
      
    if (argparser_keyword_present(pp, "-resample"))
      { o->resample = (int)argparser_get_next_int(pp, 1, nef_MAX_RESAMPLE); }
    else
      { o->resample = 1; }

    if (argparser_keyword_present(pp, "-trend"))
      { o->trend_deg = (int)argparser_get_next_int(pp, -1, nef_MAX_TREND_DEG);
        o->trend_keep = argparser_get_next_bool(pp);
      }
    else
      { o->trend_deg = -1; o->trend_keep = FALSE; }
 
    o->verbose = argparser_keyword_present(pp, "-verbose");
   
    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }
