#define PROG_NAME "nmeeg_median_filter"
#define PROG_DESC "Applies a weighted median filter to an EEG dataset."
#define PROG_VERS "2023-11-02"

#define nmeeg_median_filter_C_COPYRIGHT \
  "Copyright © 2023 by the State University of Campinas (UNICAMP)"

/* Last edited on 2023-12-07 16:51:03 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -lowWeights {NW_LO} {KIND_LO} {PARM_LO} \\\n" \
  "    -highWeights {NW_HI} {KIND_HI} {PARM_HI} \\\n" \
  "    [ -keepMean {KEEP} ] \\\n" \
  "    -outPrefix {OPREF} \\\n" \
  "    [ -writeLow ] \\\n" \
  "    [ -writeHigh ] \\\n" \
  "    " neuromat_eeg_io_FORMAT_OPT_HELP " \\\n" \
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
  " applies two running weighted median filters to it, and writes the" \
  " fitered component and/or the high- and low-frequency \"noise\" components.\n" \
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
  "  The output files will have the same format.  The {NE} output electrode" \
  " channels will be the result of the filtering applied to the" \
  " respective input electrode channels, while the non-electrode" \
  " channels will be copied without any change.\n" \
  "\n" \
  "MEDIAN FILTER\n" \
  "  The processing of each electrode channel uses two weighted median" \
  " filters, \"low-frequency\" (with a broad window of {NW_LO} samples)" \
  " and \"high-frequency\" (with a narrower window of {NW_HI} samples).\n" \
  "\n" \
  "  Each median filter takes a certain odd number {NW} of consecutive" \
  " samples, assigns to each sample a weight, and finds the value {VM} such" \
  " that the samples less than {VM} have half of the total weight.\n" \
  "\n" \
  "  The low-frequency filter, with width {NW_LO}, is" \
  " applied first, and its output {LM} is" \
  " the \"low-pass component\" or \"low-frequency residue\", presumably due" \
  " to electrode drift and jumps. This signal {LM} is" \
  " subtracted from the sample at the center of the" \
  " window, yielding an intermediate signal {TM}.  Then the" \
  " high-frequency filter, with width {NW_HI}, is applied to {TM}, yielding" \
  " the \"cleaned\" or \"band-pass component\" signal {MM}.  The" \
  " difference {TM-MM} is the \"high-pass component\" or" \
  " \"high-frequency residue\" {HM}, presumably" \
  " just meaningless high-frequency noise.\n" \
  "\n" \
  "  The window" \
  " widths are reduced near the ends of the file so that the windows" \
  " stay inside the frame index range {0..NT-1}.\n" \
  "\n" \
  "OUTPUT FILES\n" \
  "  The program always writes the cleaned signal {MM} to" \
  " file \"{OPREF}_mm.txt\".  If the options \"-writeLow\"" \
  " and/or \"-WriteHigh\" are given, the program also writes the {LM}" \
  " and/or {HM} signals, to files \"{OPREF}_lo.txt\" and" \
  " \"OPREF}_hi.txt\", respectively. \n" \
  "\n" \
  "OPTIONS\n" \
  "  -lowWeights {NW_LO} {KIND_LO} {PARM_LO}\n" \
  "  -highWeights {NW_HI} {KIND_HI} {PARM_HI}\n" \
  "    These mandatory arguments specify the width (number of samples) of the" \
  " low-frequency (broad) and high-frequency (narrow) windows" \
  " for the running median computation, respectively {NW_LO} and {NW_HI}.  Each" \
  " must be an odd number, at least 3; or zero, signifying that the" \
  " corresponding filter is not to be applied.  Namely, if {NW_LO} is zero, then" \
  " the low-frequency residue {LM} is assumed to be zero, and the intermediate" \
  " signal {TM} is the same as the input signal.  If {NW_HI} is zero, the cleaned" \
  " signal {MM} is the same as {TM}, and the high-frequency residue {HM} will be zero.\n" \
  "\n" \
  "    The {KIND_LO} and {KIND_HI} attributes are strings that specify the" \
  " general shape of the window weights. They may be:\n" \
  "" wt_table_kind_from_string_INFO "\n" \
  "    The default is \"hann\".  The arguments {PARM_LO} and {PARM_HI} are" \
  " the extra parameters required by some window kinds, shown in parentheses in the list above.\n" \
  "\n" \
  "  -keepMean {KEEP}\n" \
  "    This optional flag specifes the average (\"DC level\") of the electrodes" \
  " output files.  If {KEEP} is 1 or 'T', every electrode in every output file" \
  " will be shifted so that its mean value is the same as in the" \
  " input file.  If {KEEP} is 0 or 'F', not specified, the input mean" \
  " value will mostly go to the low-freq {LM} component, and the other two" \
  " components {MM,HM} will end up with means close to zero.  This option" \
  " does not affect marker (non-electrode) channels.  The" \
  " default is \"-keepMean F\".\n" \
  "\n" \
  "  -outPrefix {OPREF}\n" \
  "    This mandatory argument specifies the prefix for all output file names.\n" \
  "\n" \
  "  -writeLow\n" \
  "  -writeHigh\n" \
  "    These optional flags specify that the program should write" \
  " the low-frequency residue {LM} and/or the high-frequency" \
  " residue {HM}.\n" \
  "\n" \
  "" neuromat_eeg_io_FORMAT_OPT_INFO(nemf_format_DEFAULT) "\n" \
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
  "  Created 2023-11-02 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2023-11-02 Created by cannibalizing {nmeeg_filter.c}.\n" \
  "  2023-11-22 Added the high-freq filter.\n" \
  "  2023-11-23 Added the \"-keepMean\" option.\n" \
  "  2023-11-24 Added the \"-weights\" option.\n" \
  "  2023-11-28 Replaced \"width\" and \"-weights\" by \"-lowWeights\" and \"-highWeights\".\n" \
  "  2023-12-05 Added the \"-format\" option.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " nmeeg_median_filter_C_COPYRIGHT ".\n" \
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

#include <r2.h>
#include <rn.h>
#include <affirm.h>
#include <jsfile.h>
#include <argparser.h>
#include <wt_table.h>
#include <wt_table_quantize.h>
#include <wt_table_args_parse.h>
#include <wt_table_generic.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_io.h>
#include <neuromat_eeg_header.h>
#include <neuromat_eeg_source.h>
#include <neuromat_eeg_channel_stats.h>
#include <neuromat_median_filter.h>

#define nemf_MAX_WIDTH 8192
  /* Max window width for median filter. */

#define nemf_format_DEFAULT "%14.8e"
  /* Default for the "-format" option. */

typedef struct nemf_wt_options_t
  { int32_t nw;           /* Nominal window width (number of samples). */
    wt_table_kind_t kind; /* Weight distribution kind. */
    double parm;          /* Additional distribution parameter. */
  } nemf_wt_options_t;
  /* Parameters for a window weights table. */

typedef struct nemf_options_t
  { nemf_wt_options_t *lowWeights;  /* Parameters for low-freq filter. */
    nemf_wt_options_t *highWeights; /* Parameters for high-freq filter. */
    char *outPrefix;                /* Prefix for output file names. */
    bool_t writeLow;                /* True to write the low-frequency residue. */
    bool_t writeHigh;               /* True to write the high-frequency residue. */
    bool_t keepMean;                /* True to preserve the mean in all outputs. */
    char *format;                   /* Format for the data samples in output files. */
    bool_t verbose;                 /* True to print more info. */
  } nemf_options_t;
  /* Arguments from command line. */

void nemf_show_statistics(char *name, int32_t nt, int32_t nc, double **dat, int32_t ne, char *chname[]);
  /* Prints to {stderr} the per-channel and global statistics of the samples in dataset
    {dat[0..nt-1]{0..nc-1]}, assuming that the first {ne} channels are electrode potentials
    and the channel names are {chname[0..nc-1]}. */
   
int32_vec_t nemf_make_weights(nemf_wt_options_t *wop, bool_t verbose);
  /* Creates a weight table .
    
    Then converts it to integers with {wt_table_quantize} 
    with {wt_min = 0} and a large {wt_sum}. 
    
    Then eliminates zero entries at the ends of the table,
    preserving its symmetry.  Thus the length {wt.ne} of the returned vector
    may be less than {nw} elements, although it will still
    be odd and at least 3. */
    
double nemf_get_mean(int32_t ne, int32_t nt, int32_t ie, double **val);
  /* Return the average value of electrode {ie} in the array {val[0..nt-1][0..ne-1]}. */
  
void nemf_adjust_mean(int32_t ne, int32_t nt, int32_t ie, double **vot, double avg_in);
  /* Modifies the values of {vot[0..nt-1][ie]} in the array {vot[0..nt-1][0..ne-1]} 
    so that their average is {avg_in}. */

neuromat_eeg_header_t *nemf_make_header_for_output(neuromat_eeg_header_t *hin);
  /* Creates a template header suitable for the output files. The filtering steps
    shoudl be recorded with {nemf_record_median_filter_in_header} below. */

void nemf_write_eeg_signals
  ( char *pref,
    char *tag, 
    int32_t nt, 
    int32_t nc,
    double **val,
    char *fmt,
    int32_t ne,
    neuromat_eeg_header_t *hot
  );
  /* Writes the signals {val[0..nt-1][0..nc-1]}, comprising {nc}
    channels sampled at {nt} times, to file "{pref}_{tag}.txt".
    Each value is written with format {fmt}.
    
    Writes the header {hot} to the file. Checks that {hot->{nc,ne,nt}} are equal to
    {nc,ne,nt}, and assumes that these fields, including {hot->chname}, have
    been updated to reflect any channel excludions. */
  
void nemf_record_median_filter_in_header(neuromat_eeg_header_t *hot, nemf_wt_options_t *wop, bool_t med);
  /* Updates the header {hot} to record that a running median filter of width {nw}, kind {kind},
    and parameter {parm} were applied to the electrodes.  Assumes that if {med} is true
    the signal is the median value, if {med} is false the signal is the residual (original 
    minus median). */

nemf_options_t *nemf_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments.  */
  
nemf_wt_options_t *nemf_parse_wt_options(argparser_t *pp);
  /* Parses the {nw}, {kind}, and {parm} parameters 
    of a weights table.  */
  
int32_t main(int32_t argc, char **argv);
  /* Main prog. */

int32_t main(int32_t argc, char **argv)
  {
    nemf_options_t *o = nemf_parse_options(argc, argv);

    /* Read the EEG header, provide defaults: */
    int32_t nl = 0;
    neuromat_eeg_header_t *hin = neuromat_eeg_header_read(stdin, 20, 600.0, &nl);
    fprintf(stderr, "read %d header lines\n", nl);
    
    /* Get important parameters: */
    int32_t nc = hin->nc; /* Number of channels per data frame (including triggers etc.). */
    int32_t ne = hin->ne; /* Assumes that the first {ne} channels are electrode potentials. */
    fprintf(stderr, "input file has %d channels, including %d electrode potentials\n", nc, ne);

    /* Read the EEG data frames: */
    int32_t nl0 = nl;
    int32_t nt = 0;
    double **val = neuromat_eeg_data_read(stdin, 0, 0, nc, &nl, &nt);
    fprintf(stderr, "read %d lines, got %d data frames\n", nl - nl0, nt);
    demand(nt > 0, "no frames read");
    demand(nt == hin->nt, "inconsistent header {nt}");
    
    nemf_show_statistics("input signal", nt, ne, val, ne, hin->chname);

    /* Smoothed signal {vmm[0..nt-1][0..ne-1]}: */
    double **vmm = neuromat_eeg_new(nt, nc);

    /* Low and high signals, if requested: */
    double **vlo = (o->lowWeights->nw == 0 ? NULL : neuromat_eeg_new(nt, nc));
    double **vhi = (o->highWeights->nw == 0 ? NULL : neuromat_eeg_new(nt, nc));
    
    /* Copy input EEG {val} to {vmm}: */
    for (int32_t it = 0; it < nt; it++) { rn_copy(nc, val[it], vmm[it]); }
      
    /* Apply filters: */
    int32_t nw_lo, nw_hi; /* Effective filter widths. */
    if (o->lowWeights->nw != 0)
      { /* Get the low-freq component of {vmm} into {vlo}: */
        int32_vec_t wt_lo = nemf_make_weights(o->lowWeights, o->verbose);
        nw_lo = wt_lo.ne; /* Note: may be less than {o->lowWeights->nw}. */
        if (o->verbose) { fprintf(stderr, "  applying low-freq filter (nw = %d)...\n", nw_lo); }
        neuromat_median_filter_apply(nt, ne, vmm, vlo, nw_lo, wt_lo.e, o->verbose);
        for (int32_t it = 0; it < nt; it++)
          { /* Copy the non-electrode channels tp {vlo}: */
            for (int32_t ic = ne; ic < nc; ic++) { vlo[it][ic] = vmm[it][ic]; }
            /* Subtract the low-freq component from the electrode signals in {val}: */
            for (int32_t ie = 0; ie < ne; ie++) { vmm[it][ie] -= vlo[it][ie]; }
          }
       }
     else
       { nw_lo = 0; }
       
     if (o->highWeights->nw != 0)
       { /* Apply the high-pass filter on {vmm}, put the MEDIUM component temporarily in {vhi}: */
         int32_vec_t wt_hi = nemf_make_weights(o->highWeights, o->verbose);
         nw_hi = wt_hi.ne; /* Note: may be less than {o->highWeights->nw}. */
         if (o->verbose) { fprintf(stderr, "  high-freq filter (nw = %d)...\n", nw_hi); }
         neuromat_median_filter_apply(nt, ne, vmm, vhi, nw_hi, wt_hi.e, o->verbose);
         for (int32_t it = 0; it < nt; it++)
          { /* Replace {vmm,vhi} by {vhi,vmm-vhi}: */
            for (int32_t ie = 0; ie < ne; ie++) 
              { double tmm = vhi[it][ie]; 
                vhi[it][ie] = vmm[it][ie] - tmm;
                vmm[it][ie] = tmm;
              }
            /* Copy the non-electrode channels to {vhi}: */
            for (int32_t ic = ne; ic < nc; ic++) { vmm[it][ic] = val[it][ic]; vhi[it][ic] = val[it][ic]; }
          }
       } 
     else
       { nw_hi = 0; }
    
    if (o->keepMean)
      { /* Apply the mean value adjustment to the output signals: */
        if (o->verbose) { fprintf(stderr, "  adjusting means...\n"); }
        for (int32_t ie = 0; ie < ne; ie++)
          { double avg_in = nemf_get_mean(ne, nt, ie, val);
            if (o->verbose) { fprintf(stderr, "    electrode %s in = %14.8f", hin->chname[ie], avg_in); }
            nemf_adjust_mean(ne, nt, ie, vmm, avg_in);
            if (o->verbose) { fprintf(stderr, "  mm = %14.8f", nemf_get_mean(ne, nt, ie, vmm)); }
            if (nw_lo != 0)
              { nemf_adjust_mean(ne, nt, ie, vlo, avg_in);
                if (o->verbose) { fprintf(stderr, "  lo = %14.8f", nemf_get_mean(ne, nt, ie, vlo)); }
              }
            if (nw_hi != 0)
              { nemf_adjust_mean(ne, nt, ie, vhi, avg_in);
                if (o->verbose) { fprintf(stderr, "  hi = %14.8f", nemf_get_mean(ne, nt, ie, vhi)); }
              }
            if (o->verbose) { fprintf(stderr, "\n"); }
          }
       }
    
    nemf_show_statistics("cleaned signal", nt, ne, vmm, ne, hin->chname);
    if (nw_lo != 0) { nemf_show_statistics("low-freq residue", nt, ne, vlo, ne, hin->chname); }
    if (nw_hi != 0) { nemf_show_statistics("high-freq residue", nt, ne, vhi, ne, hin->chname); }

    if (TRUE)
      { /* Write the filtered output, subsampled: */
        neuromat_eeg_header_t *hot = nemf_make_header_for_output(hin);
        nemf_record_median_filter_in_header(hot, o->lowWeights, TRUE);
        nemf_record_median_filter_in_header(hot, o->highWeights, FALSE);
        nemf_write_eeg_signals(o->outPrefix, "mm", nt, nc, vmm, o->format, ne, hot); 
        neuromat_eeg_header_free(hot);
      }

    if (o->writeLow)
      { if (nw_lo == 0)
          { fprintf(stderr, "!! low-pass signal not computed, \"-writeLow\" ignored"); }
        else
          { neuromat_eeg_header_t *hot = nemf_make_header_for_output(hin);
            nemf_record_median_filter_in_header(hot, o->lowWeights, TRUE);
            nemf_write_eeg_signals(o->outPrefix, "lo", nt, nc, vlo, o->format, ne, hot); 
            neuromat_eeg_header_free(hot);
          }
      }
    
    if (o->writeHigh)
      { if (nw_hi == 0)
          { fprintf(stderr, "!! high-pass signal not computed, \"-writeHigh\" ignored"); }
        else
          { neuromat_eeg_header_t *hot = nemf_make_header_for_output(hin);
            nemf_record_median_filter_in_header(hot, o->lowWeights, FALSE);
            nemf_record_median_filter_in_header(hot, o->highWeights, TRUE);
            nemf_write_eeg_signals(o->outPrefix, "hi", nt, nc, vhi, o->format, ne, hot);
            neuromat_eeg_header_free(hot);
          }
      }
      
    neuromat_eeg_free(val, nt, nc);
    neuromat_eeg_free(vmm, nt, nc);
    if (vlo != NULL) { neuromat_eeg_free(vlo, nt, nc); }
    if (vhi != NULL) { neuromat_eeg_free(vhi, nt, nc); }
    
    return 0;
  }
      
int32_vec_t nemf_make_weights(nemf_wt_options_t *wop, bool_t verbose)
  {
    int32_t nw = wop->nw;
    demand(((nw & 1) == 1) && (nw >= 3), "invalid window size {nw}");
    int32_t hw = (nw-1)/2;
    
    /* Compute table, un-normalized: */
    wt_table_kind_t kind = wop->kind;
    double parm = wop->parm;
    double wd[nw];
    wt_table_fill(kind, nw, parm, wd, FALSE, NULL);

    /* The window weights are {wt[0..nw-1]}, where {hw_cur} is in {0..hw}: */
    int32_vec_t wt = int32_vec_new(nw);   /* Allocated with max size, but uses only {wt[0..nw_cur-1]}. */
    
    int32_t wt_min = 0;
    int32_t wt_sum = (1 << 23);
    int32_t wsum = wt_table_quantize(nw, wd, wt_min, wt_sum, wt.e);
    assert(wsum > 0);
    
    if (wt.e[0] == 0)
      { /* Eliminate zero entries: */
        int32_t hw_r = hw;
        while ((hw_r > 0) && (wt.e[hw - hw_r] == 0)) 
          { assert(wt.e[hw+hw_r] == 0);  /* Table should be symmetric. */
            hw_r--;
          }
        demand(hw_r > 0, "window weight function is too narrow");
        for (int32_t k = -hw_r; k <= +hw_r; k++)
          { wt.e[hw_r + k] = wt.e[hw + k]; }
        nw = 2*hw_r + 1;
        int32_vec_trim(&wt, nw);
        fprintf(stderr, "window size reduced to %d\n", nw);
     }
   if (verbose)
     { double wd[nw];
       for (int32_t iw = 0; iw < nw; iw++) { wd[iw] = (double)wt.e[iw]; }
       double avg = wt_table_avg(nw, wd);
       double dev = sqrt(wt_table_var(nw, wd, avg));
       fprintf(stderr, "  weight table average = %.2f deviation = %.4f\n", avg, dev);
     }
   return wt;   
  }  

double nemf_get_mean(int32_t ne, int32_t nt, int32_t ie, double **val)
  { demand((ie >= 0) && (ie < ne), "bad channel index");
    double sum = 0;
    for (int32_t it = 0; it < nt; it++)
      { double vi = val[it][ie];
        demand(isfinite(vi), "infinite or {NAN} sample");
        sum += vi;
      }
    return sum/nt;
  }

void nemf_adjust_mean(int32_t ne, int32_t nt, int32_t ie, double **vot, double avg_in)
  { demand((ie >= 0) && (ie < ne), "bad channel index");
    double avg_ot = nemf_get_mean(ne, nt, ie, vot);
    for (int32_t it = 0; it < nt; it++)
      { double vi = vot[it][ie];
        demand(isfinite(vi), "infinite or {NAN} sample");
        vot[it][ie] = vi + avg_in - avg_ot;
      }
  }

neuromat_eeg_header_t *nemf_make_header_for_output(neuromat_eeg_header_t *hin)
  { neuromat_eeg_header_t *hot = neuromat_eeg_header_new();
    neuromat_eeg_header_merge(hot, hin);
    hot->orig = neuromat_eeg_source_new();
    neuromat_eeg_header_merge_orig(hot->orig, hin->orig);
    return hot;
  }

void nemf_write_eeg_signals
  ( char *pref,
    char *tag, 
    int32_t nt, 
    int32_t nc, 
    double **val, 
    char *fmt,
    int32_t ne,
    neuromat_eeg_header_t *hot
  )
  {
    /* Set the dataset size in header: */
    hot->nt = nt;
    hot->nc = nc;
    hot->ne = ne;
    
    /* Write to {stdout}: */
    char *fname = NULL;
    asprintf(&fname, "%s_%s.txt", pref, tag);
    FILE *wr = open_write(fname, TRUE);
    neuromat_eeg_header_write(wr, hot);
    neuromat_eeg_data_write(wr, nt, nc, val, fmt, 0, nt-1, 1);
    fflush(wr);
    free(fname);
  }
    
void nemf_show_statistics(char *name, int32_t nt, int32_t nc, double **dat, int32_t ne, char *chname[])
  { neuromat_eeg_channel_stats_t *st = neuromat_eeg_channel_stats_new(nc);
    neuromat_eeg_channel_stats_t *stg = neuromat_eeg_channel_stats_new(1);
    double eps = 0.01; /* Assumed uncertainty of measurement (µV). */
    neuromat_eeg_channel_stats_gather_all(nt, nc, dat, NULL, eps, st, ne, stg);
    fprintf(stderr, "  --- channel statistics of %s ---\n", name);
    neuromat_eeg_channel_stats_print_all(stderr, 2, nc, chname, FALSE, st, ne, stg);
    fprintf(stderr, "\n");
    free(st);
    free(stg);
  }

void nemf_record_median_filter_in_header(neuromat_eeg_header_t *hot, nemf_wt_options_t *wop, bool_t med)
  {
    fprintf(stderr, "!! %s: NOT IMPLEMENTED\n", __FUNCTION__);
  }

nemf_options_t *nemf_parse_options(int32_t argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    nemf_options_t *o = talloc(1, nemf_options_t); 
    
    /* Parse keyword parameters: */
    
    argparser_get_keyword(pp, "-lowWeights");
    o->lowWeights = nemf_parse_wt_options(pp);
    
    argparser_get_keyword(pp, "-highWeights");
    o->highWeights = nemf_parse_wt_options(pp);
  
    if (argparser_keyword_present(pp, "-keepMean"))
      { o->keepMean = argparser_get_next_bool(pp); }
    else
      { o->keepMean = FALSE; }
   
    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next_non_keyword(pp);
    
    o->writeLow = argparser_keyword_present(pp, "-writeLow");
    o->writeHigh = argparser_keyword_present(pp, "-writeHigh");
  
    if (argparser_keyword_present(pp, "-format"))
      { o->format = argparser_get_next_non_keyword(pp); }
    else
      { o->format = nemf_format_DEFAULT; }
 
    o->verbose = argparser_keyword_present(pp, "-verbose");
   
    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }

nemf_wt_options_t *nemf_parse_wt_options(argparser_t *pp)
  { 
    nemf_wt_options_t *wop = talloc(1, nemf_wt_options_t); 
    
    wop->nw = (int32_t)argparser_get_next_int(pp, 0, nemf_MAX_WIDTH);
    if ((wop->nw != 0) && ((wop->nw < 3) || ((wop->nw % 2) != 1))) { argparser_error(pp, "invalid window width"); }
    
    wop->kind = wt_table_args_parse_kind(pp);
    wop->parm = argparser_get_next_double(pp, -INF, +INF); 
    
    return wop;
  }
