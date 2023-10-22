#define PROG_NAME "nmeeg_filter"
#define PROG_DESC "Applies a bandpass filter to an EEG dataset."
#define PROG_VERS "2013-11-15"

#define nmeeg_filter_C_COPYRIGHT \
  "Copyright © 2013 by the State University of Campinas (UNICAMP)"

/* Last edited on 2023-10-21 21:51:18 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    [ -rebase [ {CHAN_WT}.. ] ] \\\n" \
  "    [ -exclude {EX_NAME} ].. \\\n" \
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
  " result (possibly sub-sampled) to standard output.  Also optionally changes the potential reference to be the" \
  " average of selected channels.\n" \
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
  " of frames may be reduced and the header will include information about the filter band.  Optionally," \
  " some channels may be excluded too.\n" \
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
  "  -exclude {EX_NAME}\n" \
  "    This optional argument requests the exclusion of the channel called {EX_NAME} from" \
  " the dataset, before any further processing.\n" \
  "\n" \
  "  -filter {F_TYPE} {F_LO0} {F_LO1} {F_HI1} {F_HI0}\n" \
  "    This optional argument specifies the frequency bandpass/bandkill filter type {F_TYPE} and" \
  " the band parameters {F_LO0} {F_LO1} {F_HI1} {F_HI0}.  If omitted, no frequency filtering will be applied.\n" \
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
  "  -rebase [ {CHAN_WT}.. ] \n" \
  "    This optional argument specifies that the electrode potentials in each frame" \
  " should be replaced by the difference between the original samples and a" \
  " weighted average of all the electrode samples in the frame.  The \"-rebase\" keyword may be followed" \
  " optionally by a list of {NE} numbers {CHAN_WT}," \
  " the weights to be used for each electrode.\n" \
  "    This" \
  " correction is applied AFTER channel exclusion, trend removal, filtering, and trend" \
  " replacement. Howevr note that the weight list, if given, must include dummy weights for" \
  " any electrodes that will be excluded, even though those weights are not used.   Excluded" \
  " channels and channels with zero weight will be" \
  " ignored in the average and will not be modified. (To omit a channel from" \
  " the average but still get it rebased, use a very small weight, like 1e-10.)\n" \
  "    If the \"-rebase\" keyword is present" \
  " but is not followed by any weight, the program will assume equal" \
  " weights (simple average). If the keyword is omitted, rebasing will not be done.\n" \
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
  " discarded (\"F\" or 0).  If {TREND_DEG} is {-1} (or if the option is not" \
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
  "  2021-08-23 Added \"-rebase\".\n" \
  "  2021-08-26 Added \"-electrodes\".\n" \
  "  2021-08-30 Added \"-exclude\", removed \"-electrodes\".\n" \
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
#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <math.h>

#include <fftw3.h>

#include <r2.h>
#include <rn.h>
#include <affirm.h>
#include <jsfile.h>
#include <argparser.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_io.h>
#include <neuromat_eeg_header.h>
#include <neuromat_filter.h>
#include <neuromat_eeg_channel_stats.h>

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
  { bool_t rebase_op;          /* True iff rebasing was requested. */
    double_vec_t rebase_wt;    /* Channel weights for rebasing, or empty vector if uniform or no rebasing. */
    string_vec_t exclude;      /* Names of channels to exclude. */
    int32_t trend_deg;         /* Degree of polynomial trend to separate, or {-1} for no trend. */
    bool_t trend_keep;         /* Tells whether trend polynomial is to be kept (true) of discarded (false). */
    nef_freq_filter_t filter;  /* Freq filter parameters. */
    bool_t invert;             /* Freq filter mode: false for bandpass, true for bandkill. */
    int32_t resample;          /* Resampling step. */
    bool_t verbose;            /* True to print more info. */
  } nef_options_t;
  /* Arguments from command line. */
     
void nef_exclude_channels(int32_t nx, char *exnames[], int32_t nt, int32_t *ncP, double **val, int32_t *neP, char *chname[], double wt[]);
  /* On input, assumes that the numner of channels {nc} is {*ncP}, the number of electrodes {ne} is {*neP},
    and that {val[0..nt-1][0..nc-1]} are the sample values of the dataset. 
    Excludes from the {val} the channels whose names are in {exnames[0..nx-1]}.  Also removes them from 
    the channels name list {chname[0..nc-1]} and the rebasing weight list {wt[0..ne-1]} (if not {NULL}).
    Updates {*ncP} and {*neP} accordingly. */

void nef_rebase_frames(int32_t nt, int32_t nc, double **val, int32_t ne, double wt[]);
  /* Subtracts from every sample in the dataset {val} the average
    of all samples in the frame.  Howevr, channels with have weight 
    0 are not modified. */

void nef_filter_signals
  ( int32_t nt, 
    int32_t nc, 
    double **val, 
    int32_t ne,
    double fsmp,
    nef_freq_filter_t *ff,
    bool_t invert, 
    int32_t trdeg,
    bool_t trkeep,
    int32_t resample,
    bool_t verbose
  );
  /* Applies to each channel in the EEG dataset {val[0..nt-1][0..nc-1]} the ]
    filtering specified by {ff}, or its complement if {invert} is 1.
    Assumes that any channels to be excluded have been removed already from {ne} and {val}.

    If {trdeg} is non-negative, applies trend removal of the given degree before filtering.
    If {trkeep} is true, puts back the removed trend after filtering.
    
    The {resample} argument should be the requested resampling step. It
    is used only to check whether the requested step is consistent with
    the filtering. The dataset is NOT resampled. */
  
void nef_write_eeg_signals
  ( int32_t nt, 
    int32_t nc,
    double **val, 
    int32_t ne,
    double fsmp,
    int32_t it_ini,
    int32_t it_fin,
    int32_t it_step,
    neuromat_eeg_header_t *h
  );
  /* Writes the signals {val[0..nt-1][0..nc-1]}, comprising {nc}
    channels sampled at {nt} times, to standard output.
    Writes only one very {it_step} frames, beginning with frame
    {val[it_ini]} and ending with frame {val[it_fin]}.
    
    Writes the header {h} to the file, setting {h->nt} to the number
    {nw} of subsampled frames. Checks that {h->nc,h->ne} are equal to
    {nc,ne}, and assumes that these fields, including {h->chname}, have
    been updated to reflect any channel excludions.
    
    Also sets {h->kfmax} to a null value. Assumes that
    {fsmp} refers to the full dataset {val[0..nt-1]}, so sets {h->fsmp}
    to {fsmp/it_step}. Assumes that {h->orig->it_ini} and {h->orig->it_fin}
    describe the position of {val[0..nt-1]} in the original file, so
    adjusts those fields to account for the clipping and subsampling
    defined by {it_ini,it_fin}. */
      
void nef_update_rebase_in_header(neuromat_eeg_header_t *h, double wt[]);
  /* If {wc.ne > 0}, it should be the number of channels {N}, and {h->rebase_wt}
    should be {NULL}; then it sets {h->rebase_wt} to a vector with weights {wt[0..NC-1]}.
    If {wc.ne == 0}, does not change {h->rebase_wt}. */
  
void nef_update_trend_filter_in_header(neuromat_eeg_header_t *h, int32_t tdeg, bool_t tkeep);
  /* Updates the trend filter parameters in the header {h} to account for the
    trend processing defined by the parameters {tdeg, tkeep}. */

void nef_update_freq_filter_in_header
  ( neuromat_eeg_header_t *h, 
    double flo0, 
    double flo1, 
    double fhi1,
    double fhi0,
    int32_t finvert
  );
  /* Updates the freq filter parameters in the header {h} to account for the
    soft bandpass/bandkill filtering with parameters {flo0,flo1,fhi1,fhi0,finvert}.
    Expects the input data to be unfiltered. If the header already has
    filter parameters set, tries to merge them, but the result may not
    correctly describe the effect of the two sucessive filtering steps. */

nef_options_t *nef_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments. */
  
int32_t main(int32_t argc, char **argv);
  /* Main prog. */

#define nef_MAX_RESAMPLE 10
  /* Max downsampling step allowed (paranoia). */

#define nef_MAX_TREND_DEG 5
  /* Max downsampling step allowed (paranoia). */

int32_t main(int32_t argc, char **argv)
  {
    nef_options_t *o = nef_parse_options(argc, argv);
    
    nef_freq_filter_t *ff = &(o->filter);

    /* Read the EEG header, provide defaults: */
    int32_t nl = 0;
    neuromat_eeg_header_t *h = neuromat_eeg_header_read(stdin, 20, 600.0, &nl);
    fprintf(stderr, "read %d header lines\n", nl);
    
    /* Get important parameters: */
    int32_t nc = h->nc; /* Number of channels per data frame (including triggers etc.). */
    int32_t ne = h->ne; /* Assumes that the first {ne} channels are electrode potentials. */
    double fsmp = h->fsmp;
    fprintf(stderr, "input file has %d channels, including %d electrode potentials\n", nc, ne);
    fprintf(stderr, "input sampling frequency %.10g Hz\n", fsmp);

    /* Rebasing: */
    if (o->rebase_op)
      { if (o->rebase_wt.ne == 0)
          { /* Provide uniform weights: */
            double_vec_expand(&(o->rebase_wt), ne);
            for (int32_t ie = 0; ie < ne; ie++) { o->rebase_wt.e[ie] = 1.0; }
            double_vec_trim(&(o->rebase_wt), ne);
          }
        fprintf(stderr, "rebasing sample values with %d weights", o->rebase_wt.ne);
        for (int32_t ie = 0; ie < o->rebase_wt.ne; ie++) { fprintf(stderr, " %.4f", o->rebase_wt.e[ie]); }
        fprintf(stderr, "\n");
        demand(o->rebase_wt.ne == ne, "wrong number of rebasing weights");
      }
    
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

    /* Read the EEG data frames: */
    int32_t nl0 = nl;
    int32_t nt = 0;
    double **val = neuromat_eeg_data_read(stdin, 0, 0, nc, &nl, &nt);
    fprintf(stderr, "read %d lines, got %d data frames\n", nl - nl0, nt);
    demand(nt > 0, "no frames read");
    demand(nt == h->nt, "inconsistent header {nt}");
    
    int32_t nx = o->exclude.ne;
    if (nx > 0)
      { double *wt = (o->rebase_op ? o->rebase_wt.e: NULL);
        nef_exclude_channels(nx, o->exclude.e, nt, &nc, val, &ne, h->chname, wt);
        fprintf(stderr, "input file now has %d channels, including %d electrode potentials\n", nc, ne);
        /* Update {nc,ne} in header: */
        h->ne = ne; h->nc = nc;
        /* Update {ne} in weights: */
        if (o->rebase_op) { double_vec_trim(&(o->rebase_wt), ne); }
      }

    /* Filter if so requested: */
    if (ff->type != NULL)
      { nef_filter_signals(nt, nc, val, ne, fsmp, ff, o->invert, o->trend_deg, o->trend_keep, o->resample, o->verbose); }

    /* Rebase if so requested: */        
    if (o->rebase_op)
      { demand(o->rebase_wt.ne == ne, "invalid number of rebase weights");
        nef_rebase_frames(nt, nc, val, ne, o->rebase_wt.e);
      }
    
    /* Remember the filtering/rebasing parameters in header: */
    nef_update_trend_filter_in_header(h, o->trend_deg, o->trend_keep);
    nef_update_freq_filter_in_header(h, ff->flo0, ff->flo1, ff->fhi1, ff->fhi0, o->invert);
    if (o->rebase_op)  { nef_update_rebase_in_header(h, o->rebase_wt.e); }

    /* Subsampling parameters: */
    int32_t it_step = o->resample;
    int32_t it_ini = it_step/2;
    int32_t nw = (nt + it_step - 1 - it_ini)/it_step; /* Number of frames to write. */
    int32_t it_fin = it_ini + it_step*((nt - it_ini - 1)/it_step);
    assert((it_fin >= nt - it_step) && (it_fin < nt));
    fprintf(stderr, "saving %d frames (%d to %d step %d)\n", nw, it_ini, it_fin, it_step);

    /* Compute basic statistics of filtered electrode signals: */
    neuromat_eeg_channel_stats_t *st = neuromat_eeg_channel_stats_new(nc);
    neuromat_eeg_channel_stats_t *stg = neuromat_eeg_channel_stats_new(1);
    double eps = 0.01; /* Assumed uncertainty of measurement (µV). */
    neuromat_eeg_channel_stats_gather_all(nt, nc, val, NULL, eps, st, ne, stg);
    fprintf(stderr, "  --- channel statistics ---\n");
    neuromat_eeg_channel_stats_print_all(stderr, 2, nc, h->chname, FALSE, st, ne, stg);
    fprintf(stderr, "\n");
    free(st);
    free(stg);

    /* Write them out, subsampled: */
    nef_write_eeg_signals(nt, nc, val, ne, fsmp, it_ini, it_fin, it_step, h);
      
    int32_t it;
    for (it = 0; it < nt; it++) { free(val[it]); } 
    free(val);
    
    return 0;
  }
    
void nef_exclude_channels(int32_t nx, char *exnames[], int32_t nt, int32_t *ncP, double **val, int32_t *neP, char *chname[], double wt[])
  { int32_t nc = (*ncP);
    int32_t ne = (*neP);
    for (int32_t ix = 0; ix < nx; ix++)
      { char *name = exnames[ix];
        fprintf(stderr, "excluding channel %s...\n", name);
        int32_t ic = neuromat_eeg_find_channel_by_name(name, 0, nc-1, chname, FALSE);
        if (ic < 0) 
          { fprintf(stderr, "** invalid channel name '%s'\n", name);
            demand(FALSE, "aborted");
          }
        for (int32_t jc = ic; jc < nc-1; jc++) 
          { chname[jc] = chname[jc+1];
            for (int32_t it = 0; it < nt; it++)  { val[it][jc] = val[it][jc+1]; }
            if ((wt != NULL) && (jc < ne-1)) { wt[jc] = wt[jc+1]; }
          }
        nc--;
        if (ic < ne) { ne--; }
      }
    (*ncP) = nc;
    (*neP) = ne;
  }
  
void nef_rebase_frames(int32_t nt, int32_t nc, double **val, int32_t ne, double wt[])
  { /* Compute the total weight: */
    double sum_w = rn_sum(ne, wt);
    demand(sum_w > 1.0e-3, "weights are too small");
    
    /* Replace samples by difference to average: */
    for (int32_t it = 0; it < nt; it++)
      { double sum_wv = rn_dot(ne, wt, val[it]);
        double avt = sum_wv/sum_w;
        for (int32_t ie = 0; ie < ne; ie++) 
          { if (wt[ie] != 0) { val[it][ie] -= avt; } }
      }
  }

void nef_filter_signals
  ( int32_t nt, 
    int32_t nc, 
    double **val,
    int32_t ne, 
    double fsmp,
    nef_freq_filter_t *ff,
    bool_t invert, 
    int32_t trdeg,
    bool_t trkeep,
    int32_t resample,
    bool_t verbose
  )
  { 
    /* Filter procedures: */        
    auto complex filter_gain_unit(int32_t kf, int32_t nf, double fsmp);
    auto complex filter_gain_butterworth8(int32_t kf, int32_t nf, double fsmp);
    auto complex filter_gain_cbutterworth12(int32_t kf, int32_t nf, double fsmp);
    auto complex filter_gain_delay(int32_t kf, int32_t nf, double fsmp);
    auto complex filter_gain_gaussian(int32_t kf, int32_t nf, double fsmp);
    auto complex filter_gain_biquadratic(int32_t kf, int32_t nf, double fsmp);
    auto complex filter_gain_bisigmoid(int32_t kf, int32_t nf, double fsmp);
  
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
    if ((invert != 0) || (fcut >= fsmp/2))
      { demand(resample == 1, "no resampling allowed with bandkill filter"); }
    else 
      { int32_t max_resample = (int32_t)floor((fsmp/2)/fcut); 
        if (max_resample > nef_MAX_RESAMPLE) { max_resample = nef_MAX_RESAMPLE; }
        demand(resample <= max_resample, "resampling step too large for this filter");
      }
    
    /* Apply the filter: */
    neuromat_filter_apply(nt, ne, val, fsmp, trdeg, trkeep, filter_gain, verbose);

    return;
    
    complex filter_gain_unit(int32_t kf, int32_t nf, double fsmp)
      { 
        double w = 1.0;
        assert(ff->flo0 == 0);
        assert(ff->flo1 == 0);
        assert(ff->fhi1 == +INF);
        assert(ff->fhi0 == +INF);
        return (invert == 0 ? w : 1-w);
      }
    
    complex filter_gain_bisigmoid(int32_t kf, int32_t nf, double fsmp)
      { 
        demand((kf >= 0) && (kf < nf), "bad {kf}");
        if (2*kf > nf) { kf = nf - kf; } /* Reduce {kf} to absolute freq in {0..nf/2} */
        double f = kf*fsmp/nf;
        double wlo = neuromat_filter_lowpass_sigmoid(f, ff->flo0, ff->flo1);
        double whi = neuromat_filter_lowpass_sigmoid(f, ff->fhi1, ff->fhi0);
        double w = (1-wlo)*whi;
        return (invert == 0 ? w : 1-w);
      }
    
    complex filter_gain_butterworth8(int32_t kf, int32_t nf, double fsmp)
      { 
        double btfhi = ff->fhi1/0.8; btfhi = floor(btfhi + 0.5);   /* Highpass shoulder. */
        double btflo = ff->flo1*0.8; btflo = floor(btflo*10.0 + 0.5)/10.0; /* Lowpass shoulder. */
        int32_t btorder = 8;
        demand((kf >= 0) && (kf < nf), "bad {kf}");
        if (2*kf > nf) { kf = nf - kf; } /* Reduce {kf} to absolute freq in {0..nf/2} */
        double f = kf*fsmp/nf;
        double wlo = neuromat_filter_lowpass_butterworth(f, btflo, btorder);
        double whi = neuromat_filter_lowpass_butterworth(f, btfhi, btorder);
        double w = (1-wlo)*whi;
        return (invert == 0 ? w : 1-w);
      }
    
    complex filter_gain_cbutterworth12(int32_t kf, int32_t nf, double fsmp)
      { 
        double btfhi = ff->fhi1/0.8; btfhi = floor(btfhi + 0.5);   /* Highpass shoulder. */
        double btflo = ff->flo1*0.8; btflo = floor(btflo*10.0 + 0.5)/10.0; /* Lowpass shoulder. */
        int32_t btorder = 12;
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
        return (invert == 0 ? w : 1-w);
      }
    
    complex filter_gain_delay(int32_t kf, int32_t nf, double fsmp)
      { 
        demand((kf >= 0) && (kf < nf), "bad {kf}");
        complex w;
        if (kf == 0)
          { w = 1; }
        else if (2*kf == nf)
          { w = 0; }
        else 
          { if (2*kf > nf) { kf = kf - nf; }
            int32_t kdelay = (int32_t)floor(fsmp);
            if (2*kdelay >= nf) { kdelay = nf/4; }
            w = cexp(-kdelay*kf*2*M_PI/nf*I);
          }
        return (invert == 0 ? w : 1-w);
      }

    complex filter_gain_gaussian(int32_t kf, int32_t nf, double fsmp)
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
        return (invert == 0 ? w : 1-w);
      }        

    complex filter_gain_biquadratic(int32_t kf, int32_t nf, double fsmp)
      { 
        demand((kf >= 0) && (kf < nf), "bad {kf}");
        double w;
        if (kf == 0) 
          { w = 1; }
        else
          { double f = kf*fsmp/nf;
            w = neuromat_filter_lowpass_biquadratic(f, fsmp/2);
          }
        return (invert == 0 ? w : 1-w);
      }        
  }

void nef_write_eeg_signals
  ( int32_t nt, 
    int32_t nc, 
    double **val, 
    int32_t ne,
    double fsmp,
    int32_t it_ini,
    int32_t it_fin,
    int32_t it_step,
    neuromat_eeg_header_t *h
  )
  {
    demand((it_fin - it_ini) % it_step == 0, "invalid subsampling range/step");
    int32_t nw = (it_fin - it_ini)/it_step + 1;
    
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
    h->kfmax = INT32_MIN; /* Not meaningful for a time-domain dataset. */

    /* Write to {stdout}: */
    FILE *wr = stdout;
    neuromat_eeg_header_write(wr, h);
    neuromat_eeg_data_write(wr, nt, nc, val, it_ini, it_fin, it_step);
    fflush(wr);
  }

void nef_update_rebase_in_header(neuromat_eeg_header_t *h, double wt[])
  { assert(wt != NULL);
    demand((h->rebase_wt == NULL), "rebasing previously rebased data");
    int32_t ne = h->ne;
    h->rebase_wt = rn_alloc(ne);
    for (int32_t ie = 0; ie < ne; ie++) { h->rebase_wt[ie] = wt[ie]; }
  }

void nef_update_trend_filter_in_header(neuromat_eeg_header_t *h, int32_t tdeg, bool_t tkeep)
  {
    if (h->tdeg < 0)
      { /* Input was not trend-filtered: */
        demand(h->tkeep == INT32_MIN, "inconsistent trend filter params in header");
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
    int32_t finvert
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

nef_options_t *nef_parse_options(int32_t argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    nef_options_t *o = notnull(malloc(sizeof(nef_options_t)), "no mem"); 
    
    /* Parse keyword parameters: */
    
    o->exclude = string_vec_new(20);
    int32_t nx = 0; /* Count of channels to exclude. */
    while (argparser_keyword_present(pp, "-exclude"))
      { string_vec_expand(&(o->exclude), nx);
        o->exclude.e[nx] = argparser_get_next_non_keyword(pp);
        nx++;
      }
    string_vec_trim(&(o->exclude), nx);

    if (argparser_keyword_present(pp, "-rebase"))
      { o->rebase_op = TRUE;
        o->rebase_wt = double_vec_new(300);
        int32_t nw = 0;
        while (argparser_next_is_number(pp))
          { double_vec_expand(&(o->rebase_wt), nw);
            double wt = argparser_get_next_double(pp, 0.0000, 10.000);
            o->rebase_wt.e[nw] = wt; nw++;
          }
        double_vec_trim(&(o->rebase_wt), nw);
      }
    else
      { o->rebase_op = FALSE;
        o->rebase_wt = double_vec_new(0);
      }

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
      { o->resample = (int32_t)argparser_get_next_int(pp, 1, nef_MAX_RESAMPLE); }
    else
      { o->resample = 1; }

    if (argparser_keyword_present(pp, "-trend"))
      { o->trend_deg = (int32_t)argparser_get_next_int(pp, -1, nef_MAX_TREND_DEG);
        o->trend_keep = argparser_get_next_bool(pp);
      }
    else
      { o->trend_deg = -1; o->trend_keep = FALSE; }
 
    o->verbose = argparser_keyword_present(pp, "-verbose");
   
    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }
