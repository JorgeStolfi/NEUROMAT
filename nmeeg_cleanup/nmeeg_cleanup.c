#define PROG_NAME "nmeeg_cleanup"
#define PROG_DESC "Marks out blinks and clears bad channels from a set of runs."
#define PROG_VERS "2017-10-03"

#define nmeeg_cleanup_C_COPYRIGHT \
  "Copyright © 2017 by the State University of Campinas (UNICAMP)"
/* Last edited on 2023-11-04 01:26:17 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -inPrefix {IN_PREF} \\\n" \
  "    -runs {RUNID}.. \\\n" \
  "    -blinkTerms {NBL} \\\n" \
  "    -blinkThresh {BL_THR} \\\n" \
  "    -blinkRad {BL_RAD} \\\n" \
  "    -changeThresh {CH_THR} \\\n" \
  "    -numCorr {CORR_NM} \\\n" \
  "    -corrThresh {CORR_THR} \\\n" \
  "    -outPrefix {OUT_PREF} \\\n" \
  "    " argparser_help_info_HELP "" \

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program reads from standard input a set of filtered EEG datasets" \
  " and writes to standard output new /cleaned/ versions of the same.\n" \
  "\n" \
  "  The cleanup consists in (1) identifying a global set {GARB} of channels" \
  " that seem to be garbage in some runs, and setting them to zero; (2) identifying" \
  " blinking episodes, replacng the estimated blinking components of the" \
  " EEG by new channels \"BL0\", \"BL1\", ..., and marking off blink" \
  " episodes as \"BLMK\"; and (3) changing the potential reference to a" \
  " weighted mean of the channels and saving the shift as a new channel \"CZS\".\n" \
  "\n" \
  "   The cleanup process is iterative and starts with all frames in all" \
  " runs assumed to be blink episodes (\"BLMK\" set" \
  " to " stringify(nclup_BLMK_VALUE) ") and no channels assumed" \
  " to be garbage ({GARB} empty).  At each iteration, the program assumes" \
  " that the values of all channels in the set {GARB} are zero in all" \
  " frames, and:\n" \
  "\n" \
  "      1. performs a differential principal component analysis of all" \
  " frames in all runs, comparing currently assumed blink frames with" \
  " assumed non-blink frames, and assumes that the first {NBL} components" \
  " are potential patterns due to blinks.\n" \
  "\n" \
  "      2. removes those components from all frames, saving the coefficients" \
  " of those patterns as {NBL} new channels \"BL0\", \"BL1\", ....\n" \
  "\n" \
  "      3. identifies the frame index intervals where the Euclidean norm" \
  " of those coefficients is greater than a specified threshold.\n" \
  "\n" \
  "      4. updates \"BLMK\" to mark those intervals (suitably extended) as" \
  " blinking episodes.\n" \
  "\n" \
  "      5. for each run separately:\n" \
  "\n" \
  "          5.1. computes the covariance matrix of the {NE} electrodes in" \
  " frames NOT marked as blinks.\n" \
  "\n" \
  "          5.2. assumes that any channel that has large variance and" \
  " small correlations with other channels is bad, and adds it to a" \
  " set {GARB_RUN} private to that run, and a new global set {GARB}.\n" \
  "\n" \
                                                                                                                   "  These steps are repeated " stringify(nclup_MAX_ITER) " times or until the global set {GARB} and the \"BLMK\" channe" \
  "l seem to stabilize.\n" \
  "\n" \
 "  Then, for each frame of each run, the program computes a weighted" \
  " average {VAVG} of the electrode potentials that are not in the global {GARB} set, and subtracts that value from all channels.  The value {VAVG}, negated, is included in the frame as the new elecrode channel \"CZS\".\n" \
  "\n" \
  "INPUT FILES\n" \
  "  The input EEG data is read from files \"{IN_PREF}_r{RUNID}.txt\" where" \
  " {IN_PREF} and {RUNID} are specified in the command line.  The files should" \
  " all belong to the same subject and session.\n" \
  "\n" \
  "  The program assumes that every data line (/data frame/) in each input file contains" \
  " simultaneous samples from some fixed number {NC} of data channels, of which" \
  " the first {NE} are electrode potentials, while the rest are phase indicators (/markers/)" \
  " or other unspecified signals.  The electrode channels should not include" \
  " the reference electrode \"CZ\", which is assumed to be zero.\n" \
  "\n" \
  "OUTPUT FILES\n" \
  "  The cleaned versions of each dataset \"{IN_PREF}_r{RUNID}.txt\" is written" \
  " to files \"{OUT_PREF}_r{RUNID}.txt\".  These files will have three" \
  " additional channels:\n" \
  "\n" \
  "    \"CZS\" the reference electrode \"CZ\" (assumed to be zero) shifted like all others;\n" \
  "    \"BL0\", \"BL1\", ... the coefficients of the first {NBL} blink pattern for each frame;\n" \
  "    \"BLMK\" a marker channel that is positive in the frames where the" \
  " blink signals seem to be significant\n" \
  "\n" \
  "  The program also writes to \"{OUT_DIR}_BL0.txt\" the inferred electrode" \
  " potential pattern, as a single-frame EEG file.\n" \
  "\n" \
  "  The channels \"CZS\" and \"BL0\", \"BL1\", ... are potentials (in µV) and" \
  " thus are inserted after the first {NE} input channels.  Then the program appends" \
  " the original {NC-NE} marker channels from the input frame, unchanged.  Finally" \
  " the marker channel \"BLMK\" is" \
  " appended after those.\n" \
  "\n" \
  "  The resulting frame is then written to the output.\n" \
  "\n" \
  "  The program also writes to \"{OUT_DIR}_blinks.txt\" the intervals of time" \
  " in each run where the marker channel \"BLMK\" is positive.  The file has one" \
  " line for each run, with fields \"{SUBJ} {RUNID} {FLAG} {BL_INTS}\".  The {FLAG} currently" \
  " is always \"-\". The {BL_INTS} are zero or more pairs of times \"{T_INI}_{T_FIN}\" meaning" \
  " that the channel \"BLMK\" is positive between each {T_INI} and the corresponding {T_FIN}.  The" \
  " intervals are separated by commas.\n" \
  "\n" \
  "  The times {T_INI} and {T_FIN} are measured in seconds from the start of the run.  For" \
  " example, if channel \"BLMK\" is positive in frames 17 through 33, and the sampling" \
  " frequency is 100 Hz (one frame every 0.010 seconds), then {T_INI} will be 0.170 and {T_FIN} will" \
  " be 0.340.   In some cases, {T_INI} may be negative, and {T_FIN} may extend beyond" \
  " the end of the run.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -inPrefix {IN_PREF}\n" \
  "    This mandatory argument specifies the common prefix for the names of the input files.  It" \
  " should include the subject/session identification; for example \"flt-runs/s012\".\n" \
  "\n" \
  "  -runs {RUNID}..\n" \
  "    This mandatory argument specifies the IDs of the runs to process for this" \
  " subject/sesson, including the respective block numbers.  Each {RUNID} must be" \
  " a 5-digit decimal integer, zero padded; such as \"00305\" for run 5 in block 3.\n" \
  "\n" \
  "  -blinkTerms {NBL}..\n" \
  "    This mandatory argument specifies how many orthogonal components to assume" \
  " in the blink potential pattern.  At least one component must be specified.\n" \
  "\n" \
  "  -blinkThresh {BL_THR}..\n" \
  "    This mandatory argument specifies the threshold at which the strength of" \
  " the blink patterns is to be considered indicative of a blinking event.\n" \
  "\n" \
  "  -blinkRad {BL_RAD}\n" \
  "    Once a frame has been characterized as part of a blinking event by the" \
  " criterion \"-blinkThresh\", adjacent frames are considered part of the" \
  " blinking event too.  This mandatory argument specifies the number of frames" \
  " to incldue before and after each primary blink frame.\n" \
  "\n" \
  "  -changeThresh {CH_THR}\n" \
  "    This mandatory argument specifies the convergence criterion for the" \
  " iterative cleanup algorithm.  The cleanup is assumed to have converged" \
  " when, among other criteria, the Euclidean change in the blink pattern, between" \
  " two successive iterations, is less than {CH_THR} in all frames that have not" \
  " been included in blinking events.\n" \
  "\n" \
  "  -numCorr {CORR_NM}\n" \
  "  -corrThresh {CORR_THR}\n" \
  "    These mandatory arguments specify the criterion to detect \"garbage\" electrodes" \
  " in each run.  Namely, an electrode is considered to be \"non-garbage\" if it has" \
  " at least {COOR_NM} other electrodes with a correlation coefficient of {CORR_THR} or more.\n" \
  "\n" \
  "    This mandatory argument specifies the electrodes that are required to" \
  " qualify an electrode as \"not garbage\".\n" \
  "\n" \
  "  -outPrefix {OUT_PREF}\n" \
  "    Common prefix for output file names, including subject/session id.  For" \
  " example, \"cln-runs/s012\"." \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  neuromat_eeg_plot_signals.sh(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2017-10-03 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2017-10-03 J. Stolfi: created based on {nmeeg_comp_analysis}.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " nmeeg_cleanup_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define stringify(x) strngf(x)
#define strngf(x) #x

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <jsmath.h>
#include <affirm.h>
#include <vec.h>
#include <argparser.h>
#include <affirm.h>
#include <rmxn.h>
/* #include <rn.h> */
/* #include <jsstring.h> */
/* #include <jsfile.h> */

#include <neuromat_eeg.h>
#include <neuromat_eeg_pca.h>
#include <neuromat_eeg_channel_stats.h>

#include <nmeeg_cleanup_params.h>
#include <nmeeg_cleanup_run_data.h>
#include <nmeeg_cleanup_determine_blink_patterns.h>
#include <nmeeg_cleanup_detect_blinks.h>
#include <nmeeg_cleanup_detect_garbage_electrodes.h>
#include <nmeeg_cleanup_merge_garbage_flags.h>
#include <nmeeg_cleanup_create_output_run_data.h>

typedef struct nmeeg_cleanup_options_t
  { char *inPrefix;      /* Input file name prefix. */
    string_vec_t run;    /* IDs of blruns to process. */
    int blinkTerms;      /* Number of terms to consider in blink pattern. */
    double blinkThresh;  /* Threshold for blink/nonblink classification. */
    int blinkRad;        /* Frames to include in blinking episodes. */
    double changeThresh; /* Convergence criterion. */
    int numCorr;         /* Min number of highly correlated electrodes required. */
    double corrThresh;   /* Threshold to consider electrodes `highly correlated'. */
    char *outPrefix;     /* Ouptut file name prefix. */
  } nmeeg_cleanup_options_t;
  /* Arguments from command line. */

nmeeg_cleanup_options_t *nmeeg_cleanup_parse_options(int argc, char **argv);
  /* Parses the command line arguments. */

int main(int argc, char **argv);
  /* Main prog. */

int main(int argc, char **argv)
  {
    nmeeg_cleanup_options_t *o = nmeeg_cleanup_parse_options(argc, argv);
    bool_t debug = TRUE;
    
    int nb = o->blinkTerms;
    
    /* Read the datasets: */
    int nr = o->run.ne; /* Number of runs. */
    nmeeg_cleanup_run_data_t **run_in = notnull(malloc(nr*sizeof(nmeeg_cleanup_run_data_t *)), "no mem");
    int nc = -1;     /* Number of channels in input runs (same in all runs). */
    int ne = -1;     /* Number of electrode potentials in input runs (same in all runs). */
    int nt_max = -1; /* Max number of frames in any run. */
    for (int ir = 0; ir < nr; ir++)
      { nmeeg_cleanup_run_data_t *runi = nmeeg_cleanup_run_data_read(o->inPrefix, o->run.e[ir]);
        int nti = runi->h->nt;  /* Frames in this run. */
        int nci = runi->h->nc;  /* Channels in this run. */
        int nei = runi->h->ne;  /* Electrodes in this run. */
        /* Allocate the work fields: */
        nmeeg_cleanup_run_data_alloc_work_tables(runi, nb);
        for (int it = 0; it < nti; it++) { runi->blmk[it] = 0.0; }
        /* Record the max number of frames: */
        if (nti > nt_max) { nt_max = nti; }
        /* Save or check the number of channels {nc} and electrodes {ne}: */
        if (nc == -1) { nc = nci; } else { demand(nc == nci, "inconsistent channel counts"); }
        if (ne == -1) { ne = nei; } else { demand(ne == nei, "inconsistent electrode counts"); }
        run_in[ir] = runi;
      }
      
    /* Get the indices of the fixation and stimulus channels: */
    int nk = 2;
    int ic_mk[nk];
    char **chname_in = run_in[0]->h->chname;
    ic_mk[0] = neuromat_eeg_find_channel_by_name("FX", ne, nc-1, chname_in, TRUE);
    ic_mk[1] = neuromat_eeg_find_channel_by_name("ST", ne, nc-1, chname_in, TRUE);
    
    /* Main iteration to detect blinks and bad electrodes: */
    int iter = 0; /* Number of iterations done. */
    bool_t converged = FALSE;
    double *blP = rmxn_alloc(nb, ne);  /* Potential patterns of blinks.*/
    double *blQ = rmxn_alloc(nb, nb); /* The array {(blP*blP')^-1}. */
    bool_t garb[nc];   /* TRUE for channels that seem to be garbage in some run. */
    for (int ic = 0; ic < nc; ic++) { garb[ic] = FALSE; }
    double **work1 = alloc_C_matrix(nt_max, ne); /* Fitted blink components. */
    while ((iter < nmeeg_cleanup_MAX_ITER) && (! converged))
      { 
        fprintf(stderr, "=== BEGIN ITERATION %d =====================================\n", iter);
        
        /* Determine the blink patterns: */
        fprintf(stderr, "--- determining blink patterns --------------------------------\n");
        bool_t onlyBlinks = (iter > 0); /* Consider everybody in iteration 0. */
        nmeeg_cleanup_determine_blink_patterns(nr, nc, ne, run_in, onlyBlinks, garb, nb, blP);

        /* Compute the linear system's matrix {blQ}: */
        neuromat_eeg_pca_compute_fitting_matrix(nb, ne, blP, blQ);

        if (debug)
          { fprintf(stderr, "  the {blP} matrix:\n");
            for (int ie = 0; ie < ne; ie++)
              { fprintf(stderr, "    channel %3d = %-6s", ie, chname_in[ie]);
                for (int ib = 0; ib < nb; ib++)
                  { fprintf(stderr, " %10.4f", blP[ib*ne + ie]); }
                fprintf(stderr, "\n");
              }
            fprintf(stderr, "  the {blQ} matrix:\n");
            for (int ib = 0; ib < nb; ib++)
              { fprintf(stderr, "    ");
                for (int jb = 0; jb < nb; jb++)
                  { fprintf(stderr, " %10.4f", blQ[ib*nb + jb]); }
                fprintf(stderr, "\n");
              }
          }

        /* Detect blinks and garbage: */
        converged = TRUE;
        for (int ir = 0;  ir < nr; ir++)
          { nmeeg_cleanup_run_data_t *runi = run_in[ir];
            
            /* Detect blinks in this run: */
            fprintf(stderr, "--- detecting blinks in run %d = %s -----------------\n", ir, runi->runid);
            bool_t convi_blink = nmeeg_cleanup_detect_blinks
              ( runi, nb, blP, blQ,
                o->blinkThresh, nmeeg_cleanup_BLMK_VALUE, o->blinkRad, o->changeThresh, 
                work1
              );
            fprintf(stderr, "  %s ", runi->runid);
            nmeeg_cleanup_run_data_print_blinks(stderr, runi, FALSE);
            fprintf(stderr, " = ");
            nmeeg_cleanup_run_data_print_blinks(stderr, runi, TRUE);
            fprintf(stderr, "\n");
              
            /* Detect garbage channels in this run: */
            fprintf(stderr, "--- detecting garbage channels in run %d = %s -----------------\n", ir, runi->runid);
            bool_t skipBlinks = TRUE; /* Ignore blink episodes for garbage channel detection. */
            bool_t convi_garb = nmeeg_cleanup_detect_garbage_electrodes
              ( runi, nk, ic_mk, skipBlinks, o->numCorr, o->corrThresh );
            
            /* Check for convergence: */
            if ((! convi_blink) || (! convi_garb)) { converged = FALSE; }
          }
          
        /* Unite all garbage channels found for next iteration: */
        bool_t conv_garb = nmeeg_cleanup_merge_garbage_flags(nr, run_in, ne, garb);
        if (! conv_garb) { converged = FALSE; }
        
        fprintf(stderr, "=== END ITERATION %d =====================================\n", iter);
        iter++;
      }
      
    /* Write cleaned runs: */
    for (int ir = 0;  ir < nr; ir++)
      { nmeeg_cleanup_run_data_t *runi = run_in[ir];
        nmeeg_cleanup_run_data_t *runo = nmeeg_cleanup_create_output_run_data(runi, ne, garb);
        nmeeg_cleanup_run_data_write(o->outPrefix, runo->runid, runo);
        if (debug) { fprintf(stderr, "freeing {runo}\n"); }
        nmeeg_cleanup_run_data_free(runo);
        if (debug) { fprintf(stderr, "freeing {runi}\n"); }
        nmeeg_cleanup_run_data_free(runi);
      }
    if (debug) { fprintf(stderr, "freeing {blP,blQ}\n"); }
    free(blP); free(blQ);
    if (debug) { fprintf(stderr, "freeing {work1}\n"); }
    free_C_matrix(work1, nt_max);
    return 0;
  }

nmeeg_cleanup_options_t *nmeeg_cleanup_parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    /* Allocate the command line argument record: */
    nmeeg_cleanup_options_t *o = notnull(malloc(sizeof(nmeeg_cleanup_options_t)), "no mem");

    /* Parse keyword parameters: */
    argparser_get_keyword(pp, "-inPrefix");
    o->inPrefix = argparser_get_next_non_keyword(pp);
    
    argparser_get_keyword(pp, "-runs");
    int nr = 0; /* Number of runs.*/
    o->run = string_vec_new(200);
    while (argparser_next_is_non_keyword(pp))
      { string_vec_expand(&(o->run), nr);
        o->run.e[nr] = argparser_get_next_non_keyword(pp);
        nr++;
      }
    string_vec_trim(&(o->run), nr);

    argparser_get_keyword(pp, "-blinkTerms");
    o->blinkTerms = (int)argparser_get_next_int(pp, 1, 1000);

    argparser_get_keyword(pp, "-blinkThresh");
    o->blinkThresh = argparser_get_next_double(pp, 0.01, 1000.0);

    argparser_get_keyword(pp, "-blinkRad");
    o->blinkRad = (int)argparser_get_next_int(pp, 0, 10000);

    argparser_get_keyword(pp, "-changeThresh");
    o->changeThresh = argparser_get_next_double(pp, 0.0001, 1000.0);

    argparser_get_keyword(pp, "-numCorr");
    o->numCorr = (int)argparser_get_next_int(pp, 1, 256);

    argparser_get_keyword(pp, "-corrThresh");
    o->corrThresh = argparser_get_next_double(pp, 0.0, 1.0);

    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next_non_keyword(pp);
    
    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }

