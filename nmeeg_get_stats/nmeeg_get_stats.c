#define PROG_NAME "nmeeg_get_stats"
#define PROG_DESC "Obtains channel statistics from an EEG dataset."
#define PROG_VERS "2013-11-15"

#define nmeeg_get_stats_C_COPYRIGHT \
  "Copyright © 2021 by the State University of Campinas (UNICAMP)"

/* Last edited on 2023-12-06 07:05:37 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -outPrefix {OPREF} \\\n" \
  "    -hStep {HSTEP} \\\n" \
  "    [ -verbose ] \\\n" \
  "    " argparser_help_info_HELP " \\\n" \
  "    < {INFILE} \\\n" \
  "    > {STATS_FILE}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program reads from standard input an EEG dataset in the plaintext (\".txt\") format," \
  " computes the per-channel and global ststistics, and writes them to standard output.\n" \
  "\n" \
  "  The program also writes a file \"{OPREF}{CHNAME}_hist.txt\" for each electrode, with the" \
  " histogram of channel values; where {CHNAME} is the electrode's name, as in the input" \
  " file header.  Each line of the file has fields \"{VLO} {VHI} {COUNT}\" where {VLO,VHI} are" \
  " two electrode potential values and {COUNT} is the count of samples of that electrode" \
  " that fall in the range {[VLO_VHI]}.\n" \
  "\n" \
  "  The input file should contain a certain" \
  " number {NT} of data frames, each with the same number {NC} of channel" \
  " values, including {NE} electrode readings and {NC-NE} non-electrode" \
  " channels (such as trigger and phase markers).  The per-channel statistics will" \
  " include the markers, but the global statistics will only look at the eletrodes.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -outPrefix {OPREF}\n" \
  "    This mandatory argument specifies the common prefix for output file names.\n" \
  "\n" \
  "  -hStep {HSTEP}\n" \
  "    This mandatory argument specifies the width of each bin of the histogram.\n" \
  "\n" \
  "  -verbose\n" \
  "    This optional keywork requests extra disgnostic messages.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  neuromat_eeg_plot_signals.sh(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2021-09-02 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2021-09-02 Created by heavily trimming down {nmeeg_filter.c}.\n" \
  "  2023-11-02 Added \"-outPrefix\", \"-hStep\", and histogram output.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " nmeeg_get_stats_C_COPYRIGHT ".\n" \
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
#include <neuromat_eeg_channel_stats.h>

typedef struct ngs_options_t
  { bool_t verbose;       /* Be verbose. */
    char *outPrefix;      /* Prefix for output file names. */
    double hStep;         /* Width of histogram bins. */
  } ngs_options_t;
  /* Arguments from command line. */
     
ngs_options_t *ngs_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments. */
  
int32_t main(int32_t argc, char **argv);
  /* Main prog. */

int32_t main(int32_t argc, char **argv)
  {
    ngs_options_t *o = ngs_parse_options(argc, argv);
    if (o->verbose) {  fprintf(stderr, "starting...\n"); }
    
    /* Read the EEG header, provide defaults: */
    int32_t nl = 0;
    neuromat_eeg_header_t *h = neuromat_eeg_header_read(stdin, 20, 600.0, &nl);
    /* fprintf(stderr, "read %d header lines\n", nl); */
    
    /* Get important parameters: */
    int32_t nc = h->nc; /* Number of channels per data frame (including triggers etc.). */
    int32_t ne = h->ne; /* Assumes that the first {ne} channels are electrode potentials. */
    fprintf(stderr, "input file has %d channels, including %d electrode potentials\n", nc, ne);

    /* Read the EEG data frames: */
    int32_t nl0 = nl;
    int32_t nt = 0;
    double **val = neuromat_eeg_data_read(stdin, 0, 0, nc, &nl, &nt);
    fprintf(stderr, "read %d lines, got %d data frames\n", nl - nl0, nt);
    demand(nt > 0, "no frames read");
    demand(nt == h->nt, "inconsistent header {nt}");
    
    /* Compute statistics of electrode signals: */
    neuromat_eeg_channel_stats_t *st = neuromat_eeg_channel_stats_new(nc);
    neuromat_eeg_channel_stats_t *stg = neuromat_eeg_channel_stats_new(1);
    double eps = 0.01; /* Assumed uncertainty of measurement (µV). */
    neuromat_eeg_channel_stats_gather_all(nt, nc, val, NULL, eps, st, ne, stg);
    fprintf(stdout, "input file has %d frames and %d channels, including %d electrode potentials\n", nt, nc, ne);
    neuromat_eeg_channel_stats_print_all(stdout, 2, nc, h->chname, FALSE, st, ne, stg);
    fflush(stdout);
    
    /* Compute electrode histograms: */
    demand(h->chname != NULL, "no channel names in file header");
    for (int32_t ie = 0; ie < ne; ie++)
      { /* Determine number {nh} of bins in histogram and allocate table {hist[0..nh-1]}: */
        double avg = st[ie].avg;
        int32_t ih_min = (int32_t)floor((st[ie].min - avg)/o->hStep - 0.001);
        int32_t ih_max = (int32_t)ceil((st[ie].max - avg)/o->hStep + 0.001);
        int32_t nh = ih_max - ih_min + 1;
        int32_t *hist = talloc(nh, int32_t);
        /* Collect the histogram of channel {ie} in {hist[0..nh-1]}: */
        for (int32_t it = 0; it < nt; it++) 
          { double valet = val[it][ie];
            demand(! isnan(valet), "NaN value in electrode");
            int32_t ih = (int32_t)floor((valet - avg)/o->hStep + 0.5);
            assert((ih >= ih_min) && (ih <= ih_max));
            hist[ih - ih_min]++;
          }
        /* Write the histogram out: */
        char *hname = NULL;
        asprintf(&hname, "%s%s_hist.txt", o->outPrefix, h->chname[ie]);
        FILE *hwr = open_write(hname, o->verbose);
        int32_t hsum = 0;
        for (int32_t ih = ih_min; ih <= ih_max; ih++)
          { double vlo = avg + (ih - 0.499999)*o->hStep;
            double vhi = avg + (ih + 0.499999)*o->hStep;
            fprintf(hwr, "%+22.6f %+22.6f %7d\n", vlo, vhi, hist[ih - ih_min]);
            hsum += hist[ih - ih_min];
          }
        assert(hsum == nt);
        fclose(hwr);
        free(hname);
        free(hist);
      }
      
    for (int32_t it = 0; it < nt; it++) { free(val[it]); } 
    free(val);
    free(st);
    free(stg);
   
    if (o->verbose) {  fprintf(stderr, "done.\n"); }

    return 0;
  }
    
ngs_options_t *ngs_parse_options(int32_t argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    ngs_options_t *o = notnull(malloc(sizeof(ngs_options_t)), "no mem"); 
    
    /* Parse keyword parameters: */

    o->verbose = argparser_keyword_present(pp, "-verbose");

    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next_non_keyword(pp);

    argparser_get_keyword(pp, "-hStep");
    o->hStep = argparser_get_next_double(pp, 0.001, 1.0e+9);

    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }
