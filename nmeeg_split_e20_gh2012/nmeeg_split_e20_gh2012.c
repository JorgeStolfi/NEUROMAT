#define PROG_NAME "nmeeg_split_e20_gh2012"
#define PROG_DESC "Splits a rawtext multiple-run 20-electrode EEG data file into individual runs."
#define PROG_VERS "2013-06-05"

#define nmeeg_split_e20_gh2012_C_COPYRIGHT \
  "Copyright © 2013 by the State University of Campinas (UNICAMP)"

/* Last edited on 2023-10-21 21:49:09 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -subject {SUBJECT} \\\n" \
  "    [ -firstRun {FRUN} ] \\\n" \
  "    -channels {CHNAME[0]} {CHNAME[1]} ... {CHNAME[NC-1]} \\\n" \
  "    [ -nElectrodes {NE} ] \\\n" \
  "    -runTypes {RTYPE[0]} {RTYPE[1]} ... {RTYPE[NR-1]} \\\n" \
  "    -trigger {TRIGNAME} \\\n" \
  "    -pulsesPerRun {NP} \\\n" \
  "    -fSampling {FSMP} \\\n" \
  "    [ -framesPerRun {FPR} ] \\\n" \
  "    [ -marker  ] \\\n" \
  "    -source {SOURCE_NAME} \\\n" \
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
  "  The program reads from standard input a 20-electrode raw EEG data file" \
  " from the Bio/NonBio experiment by Ghislain et al [Gh2013a], and extracts from" \
  " it the segments that correspond to individual experimental runs, as separate files.\n" \
  "\n" \
  "INPUT FORMAT\n" \
  "  The input file should be a raw (unfiltered an uncut) EEG recording" \
  " of the 20-electrode equipment in the \".txt\" format.\n" \
  "\n" \
  "  The program assumes that the input file contains a certain" \
  " number {NT} of data frames, each with the same number {NC} of channel" \
  " readings.\n" \
  "\n" \
  "  One of the channels in the input file should contain trigger pulses that mark" \
  " specific events in each run.  The program assumes that each run" \
  " spans {NP} consecutive (but separate) pulses in this channel.  More" \
  " precisely, each run is assumed to begin" \
  " with the first data frame of its first trigger pulse and extend" \
  " up to and including the last frame of its last trigger pulse.\n" \
  "\n" \
  "OUTPUT FILES\n" \
  "  The output files are written with" \
  " names \"{OUT_DIR}/s{SUBJECT}_r{RUN}.txt\" where {SUBJECT} is" \
  " specified in the command line with the \"-subject\" option, {RUN} is" \
  " a sequential run number starting with {FRUN}.\n" \
  "\n" \
  "  Each output file contains the data frames spanned by the {NP} marker" \
  " pulses of the corresponding run, plus" \
  " some context data frames before and after them.  Optionally the" \
  " program may adjust the amount of context data extracted at" \
  " both ends of the run to complete a fixed" \
  " number {FPR} of frames for every run.\n" \
  "\n" \
  "  Each output file will contain, in addition to the orginal {NC} channels, a certain" \
  " number {NM} of additional `marker' channels" \
  " that identify specific segments of the run tied to the marker" \
  " pulses.  These extra channels are specified by the \"-marker\" options.\n" \
  "\n" \
  "  Each output file" \
  " has some header records that" \
  " specify some relevant parameters such as number of frames" \
  " and channels, sampling frequency, run type (bio or non-bio), etc..\n" \
  "\n" \
  "OPTIONS\n" \
  "  -subject {SUBJECT}\n" \
  "    This mandatory argument is is the subject's ID number.  It is saved" \
  " in the ouput file headers and also used to form output file names.\n" \
  "\n" \
  "  -firstRun {FRUN}\n" \
  "    This optional argument specifies the ID number of the first" \
  " run in the input file, among all runs of the subject.  If" \
  " omitted, the program assumes {FRUN=1}.\n" \
  "\n" \
  "  -channels {CHNAME[0]} {CHNAME[1]} ... {CHNAME[NC-1]}\n" \
  "    This mandatory argument specifies the number {NC} and the" \
  " names {CHNAME[0..NC-1]} of the the channels in the raw data file, in" \
  " the order of the corresponding values in each data frame.  Each name" \
  " must start with an ASCII letter and contain only ASCII" \
  " letters, digits 0-9, periods '.' and underscores '_'.\n" \
  "\n" \
  "  -trigger {TRIGNAME}\n" \
  "    This mandatory argument specifies the name {TRIGNAME} of the trigger channel that" \
  " marks the boundaries of runs and run phases.  The channel is expected to" \
  " be zero most of the time and positive during the boundary pulses, which" \
  " may last for one or more consecutive frames.  The {TRIGNAME} must be" \
  " one of the strings {CHNAME[0..NC-1]}.\n" \
  "\n" \
  "  -nElectrodes {NE}\n" \
  "    This optional argument specifies the number NE} of channels that" \
  " contain electrode potentials.  They are assumed to be the first {NE} channels" \
  " in each data frame ({CHNAME[0..NE-1]}) and must not include the" \
  " trigger channel.  The parameter {NE} is saved in the output file" \
  " headers for documentation.  If this argument is omitted, the program" \
  " assumes all channels before the trigger are electrode potentials.\n" \
  "\n" \
  "  -pulsesPerRun {NP}\n" \
  "    This mandatory argument specifies the number {NP} of trigger pulses in" \
  " each experimental run. A run with {NP} pulses is divided into {NP+1} phases, namely" \
  " a `preamble' before the first pulse, {NP-1} inter-pulse phases, and" \
  " a `postamble' after the last pulse.\n" \
  "\n" \
  "  -fSampling {FSMP}\n" \
  "    This mandatory argument specifies the data sampling rate in" \
  " hertz (data frames per second).  It is recorded in the output file" \
  " headers and used to convert frame counts to time spans.\n" \
  "\n" \
  "  -framesPerRun {FPR}\n" \
  "    This optional argument specifies the number of frames to extract" \
  " for each run.  It must be at least {NPRMAX+2} where {NPRMAX} is the" \
  " maximum number of frames" \
  " actually part of any run, as delimited by the trigger channel pulses.  If a" \
  " specific run actually has {NPR < FPR} frames, the program will" \
  " extract {FPR - NPR} additional `context' frames around it, so" \
  " as to leave the run centered in the extracted segment.  If this" \
  " argument is omitted, a fixed number of context frames" \
  " is extracted at each end.\n" \
  "\n" \
  "  -marker {M_NAME} {M_IP} {M_EDGE} {M_SKIP} {M_NT} {M_VAL}\n" \
  "    Each occurrence of this optional argument specifies that" \
  " another marker channel called \"{M_NAME}\" should be added" \
  " to each data frame, besides the trigger channel.  This channel" \
  " will be set to the value {M_VAL} on a {M_NT} consecutive frames" \
  " within the run, tied to pulse number {M_IP} of the run (counting" \
  " from 0 to {NP-1}). The first marked frame will be located {M_SHIFT} frames" \
  " ahead of the first frame (if {M_EDGE=0}) or last frame (if {M_EDGE=1}) of" \
  " that trigger pulse.  The name should follow the same restrictions" \
  " as existing channel names.  The arguments {M_NT} and {M_VAL} must be" \
  " positive; argument {M_SKIP} may be zero or negative.  In any case" \
  " the marked segment must lie entirely within the extracted run.\n" \
  "\n" \
  "  -runTypes {RTYPE[0]} {RTYPE[1]} ... {RTYPE[NR-1]}\n" \
  "    This optional argument specifies a `type' string for each run" \
  " expected to be in the input file, which is recorded in the corresponding" \
  " output file header.  Each type {RTYPE[i]} must start with an ASCII letter" \
  " and contain only ASCII letters, digits 0-9, periods '.' and" \
  " underscores '_'.  If this argument is omitted, the \"type\" entry" \
  " in the output file headers is omitted.\n" \
  "\n" \
  "  -source {SOURCE_NAME}\n" \
  "    This mandatory argument specifies the source file" \
  " name (e.g. \"s1_bl2\"), with or without directories and/or" \
  " extensions.  It is used only for documentation, in the output" \
  " file headers.\n" \
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
  "REFERENCES\n" \
  "  [Gh2013a] \"Electrophysiological correlates of biological motion permanence" \
  " in humans.\" Ghislain Saunier, Eduardo F. Martins, Elisa C. Dias, José M. de" \
  " Oliveira, Thierry Pozzo, Claudia D. Vargas. Behavioural Brain" \
  " Research, volume 236 (2013) 166--174\n" \
  "AUTHOR\n" \
  "  Created 2013-06-03 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2013-06-05 Converted to {argparser.h} by J. Stolfi, IC-UNICAMP.\n" \
  "  2013-06-08 Changed \"-outPrefix\" to \"-outDir\".\n" \
  "  2013-06-10 Added \"-marker\" option.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " nmeeg_split_e20_gh2012_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <math.h>

#include <argparser.h>
#include <vec.h>
#include <affirm.h>
#include <jsfile.h>
#include <jsmath.h>
#include <jsstring.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_io.h>
#include <neuromat_eeg_header.h>

typedef struct nes_marker_t
  { char *name;  /* Channel name. */
    int ip;      /* Pulse index in run. */
    int edge;    /* Pulse edge frame: 0 (first frame) or 1 (last frame). */
    int shift;   /* Signed displacement from pulse edge. */
    int nt;      /* Number of frames to mark. */
    double val;  /* Value for marked frames. */
  } nes_marker_t;
  /* Specification for an extra `run phase' marker channel. */
  
vec_typedef(nes_marker_vec_t,nes_marker_vec,nes_marker_t);

typedef struct nes_options_t
  { 
    int subject;       /* ID number of subject (for documentation). */
    int firstRun;      /* External number of first run found in input file. */
    int framesPerRun;  /* Desired number of data frames per run, or 0 if free. */
    int nElectrodes;   /* Num of channels that are electrode potentials. */
    int trigger_ic;    /* Index of trigger channel in {channels}. */
    int pulsesPerRun;  /* Number of trigger pulses per run. */
    char *source;      /* Name of raw input file (for documentation). */
    char *outDir;      /* Directory for all output files. */
    double fSampling;  /* Assumed sampling frequency (Hz). */
    nes_marker_vec_t markers;  /* Marker channel specs. */
    string_vec_t chname;      /* Names of the channels. */
    string_vec_t runTypes;     /* Type of each run, or an empty vector. */
  } nes_options_t;
  /* Arguments from command line. */

void nes_write_eeg_run
  ( char *outDir, 
    neuromat_eeg_header_t *h, 
    int nt, 
    int nc, 
    double **vru, 
    int np, 
    int it_p_ini[], 
    int it_p_fin[]
  );
  /* Writes an extracted segment of the eeg signal array
    {vru[0..nt-1][0..nc-1]}, comprising {nc} channels sampled at {nt}
    times, to the file
    "{outDir}/s{h->orig->subject}_r{h->orig->run}.txt". The segment is
    assumed to span the data frames of the original dataset {val} with indices
    {h->orig->it_ini} to {h->orig->it_fin}, inclusive.
    
    Assumes that the run contains {np} trigger pulses, and the pulse
    with index {k} (in {0..np-1}) starts with frame {it_p_ini[k]} and
    ends with frame {it_p_fin[k]} of the original file. (Currently this
    data is not written to the file, only printed on {stderr}.)
    
    Before the actual data, also writes some header lines with the
    non-null information from {h}. Then calls
    {neuromat_eeg_data_write} to write the data. */

void nes_print_run_info(FILE *wr, char *class, int subject, int ir, char *type, int it_ini, int it_fin, double fsmp);
  /* Prints to {wr} a one-line summary for a run {ir} (external
    numbering) of the specified {subject}. Assumes that the run has the
    given {type} and spans from frame {it_ini} to frame {it_fin} of the
    original data file. The string {class} is printed in the summary.
    The sampling frequency {fsmp} is used to convert frame indices into
    time coordinates. */

void nes_print_run_pulses(FILE *wr, int np, int it_p_ini[], int it_p_fin[], int it_ini, int it_fin, double fsmp);
  /* Prints to {wr} the basic info about the {np} pulses in a run.
    Assumes that the first frame of pulse {ip} (in {0..np-1}) of the run
    is frame {it_p_ini[ip]} of the original file, and the last one is
    {it_p_fin[ip]}. Assumes that the run starts with frame number
    {it_ini} and ends with frame {it_fin} of the original data file. The
    string {class} is printed in the summary. The sampling frequency
    {fsmp} is used to convert frame indices into time coordinates. */

nes_options_t *nes_parse_options(int argc, char **argv);
  /* Parses the command line arguments. */

string_vec_t nes_parse_next_token_args(argparser_t *pp);
  /* Parses the next command line arguments a list of tokens 
    such as channel names or run types. */ 
  
int main(int argc, char **argv);
  /* Main prog. */

int main(int argc, char **argv)
  {
    nes_options_t *o = nes_parse_options(argc, argv);
    
    int nc_in = o->chname.ne;      /* Number of channels in input file. */
    int np_run = o->pulsesPerRun;   /* Number of pulses per run. */
    double fsmp = o->fSampling;     /* Sampling frequency in Hz (samples per second). */
    int ict = o->trigger_ic;        /* Index of trigger channel in input data frames (from 0). */
    int ne = o->nElectrodes;        /* Channels {0..ne-1} are electrode potentials. */
    
    int ic;
    
    fprintf(stderr, "original file name %s\n", o->source);
    fprintf(stderr, "subject ID %03d\n", o->subject);
    fprintf(stderr, "sampling frequency = %.10g\n", fsmp);

    fprintf(stderr, "assuming %d input channels = ", nc_in);
    for (ic = 0; ic < nc_in; ic++) { fprintf(stderr, " %s", o->chname.e[ic]); }
    fprintf(stderr, "\n");
    fprintf(stderr, "channels 0..%d (%d channels) are electrodes\n", ict, ne);
    fprintf(stderr, "trigger channel is number %d = %s\n", ict, o->chname.e[ict]);
    fprintf(stderr, "assuming %d trigger pulses per run\n", np_run);

    fprintf(stderr, "first experimental run in input file is number %03d\n", o->firstRun);
    if (o->runTypes.ne > 0) { fprintf(stderr, "expecting %d experimental runs in input file\n", o->runTypes.ne); }
    if (o->framesPerRun > 0) { fprintf(stderr, "writing %d data frames per run\n", o->framesPerRun); }
    fprintf(stderr, "output files are named \"%s/s%03d_r{NNN}.txt\"\n", o->outDir, o->subject);
   
    /* Read the EEG data: */
    int nl = 0; /* Number of next line from {stdin}. */
    int nt = 0; /* Number of data frames read. */
    double **val = neuromat_eeg_data_read(stdin, 0, 0, nc_in, &nl, &nt);
    fprintf(stderr, "read %d frames from %d lines\n", nt, nl-1);
    
    /* Create the names of output channels, as newly allocated copies: */
    int nm = o->markers.ne;
    int nc_out = nc_in + nm; /* Original channels plus new marker channels. */
    char **chname_out = notnull(malloc(nc_out*sizeof(char*)), "no mem"); 
    for (ic = 0; ic < nc_in; ic++) { chname_out[ic] = txtcat("", o->chname.e[ic]); }
    /* Append the names of the new marker channels: */
    int im;
    for (im = 0, ic = nc_in; im < nm; im++, ic++) { chname_out[ic] = txtcat("", o->markers.e[im].name); }
    assert(ic == nc_out);

    /* Estimate the max number of pulses {np_max}: */
    int nt_pulse_min = (int)ceil(0.1/fsmp); /* Assumes pulses repeat at most 10 times per second. */
    int np_max = nt/nt_pulse_min + 1; /* Max pulses expected in file. */
    
    /* Count and locate the pulses: */
    demand((0 <= ict) && (ict < nc_in), "invalid trigger index");
    double vmin_trig = 0;    /* "Low" trigger channel value. */
    double vmax_trig = 1000; /* "High" trigger channel value. */
    int *it_pulse_ini = notnull(malloc(np_max*sizeof(int)), "no mem"); /* Indices of trigger-up frames in {0..nt-1}. */
    int *it_pulse_fin = notnull(malloc(np_max*sizeof(int)), "no mem"); /* Indices of trigger-down frames in {0..nt-1}. */
    int np_tot = neuromat_eeg_locate_pulses(nt, nc_in, val, ict, vmin_trig, vmax_trig, np_max, it_pulse_ini, it_pulse_fin);
    demand(np_tot > 0, "no pulses detected in input file");
    
    /* Compute the number {nr} of runs in the dataset: */
    int nr = np_tot/np_run; /* Number of runs in file. */
    fprintf(stderr, "detected %d pulses (%d runs) in input file\n", np_tot, nr);
    demand(np_tot % np_run == 0, "invalid trigger pulse count");
      
    /* Decide the ideal number of frames to take before and after any run (if {o->framesPerRun == 0}): */
    double ideal_ctx_pre =  1.25; /* Ideal context stretch (seconds)  before each run. */
    double ideal_ctx_pos =  1.25; /* Ideal context stretch (seconds) after each run. */ 
    int nt_ctx_pre = (int)floor(fsmp*ideal_ctx_pre + 0.5); /* Ideal context frames before run. */
    int nt_ctx_pos = (int)floor(fsmp*ideal_ctx_pos + 0.5); /* Ideal context frames after run. */
    if (o->framesPerRun == 0) 
      { fprintf(stderr, "context frames taken before each run = %d (%.4fs)\n", nt_ctx_pre, nt_ctx_pre/fsmp);
        fprintf(stderr, "context frames taken after each run =  %d (%.4fs)\n", nt_ctx_pos, nt_ctx_pos/fsmp);
      }
    
    /* Create a header record with fields common to all runs: */
    neuromat_eeg_header_t *h = neuromat_eeg_header_new(); /* Header for run files. */
    h->orig->file = txtcat(o->source, "");
    h->orig->nt = nt;
    h->orig->fsmp = fsmp;
    h->orig->subject = o->subject;
    
    h->nc = nc_out;
    h->ne = ne;
    h->chname = chname_out;
    h->fsmp = fsmp;
    
    /* Write expanded runs: */
    int ir;
    for (ir = 0; ir < nr; ir++)
      { 
        /* Get the tight frame index range of run: */
        int it_ini = it_pulse_ini[ir*np_run]; /* Index of first frame of run in {0..nt-1}. */
        int it_fin = it_pulse_fin[(ir+1)*np_run - 1]; /* Index of last frame of run in {0..nt-1}. */
        char *run_type = o->runTypes.e[ir];
        nes_print_run_info(stderr, "strict", o->subject, o->firstRun + ir, run_type, it_ini, it_fin, fsmp);
       
        /* Adjust the frame index range {it_ini..it_fin}: */
        if (o->framesPerRun == 0) 
          { /* Run length is flexible; add fixed context frames before and after: */
            it_ini = it_ini - nt_ctx_pre;
            it_fin = it_fin + nt_ctx_pos;
          }
        else
          { /* Run length is fixed; divide the extra frames alloted between the two ends: */
            int nt_ctx = o->framesPerRun - (it_fin - it_ini + 1);
            demand(nt_ctx > 0, "cannot honor {o->framesPerRun}, actual run is too long");
            it_fin = it_fin + nt_ctx/2;
            it_ini = it_fin - o->framesPerRun + 1;
          }

        /* Check whether the context frames exist and are indeed idle: */
        demand(it_fin < (ir == nr-1 ? nt : it_pulse_ini[(ir+1)*np_run]), "not enough idle frames after end of run");
        demand(it_ini >= (ir == 0 ? 0 : it_pulse_fin[ir*np_run-1]+1), "not enough idle frames before start of run");

        /* Frame count in this run: */
        int nt_run = it_fin - it_ini + 1;
        nes_print_run_info(stderr, "padded", o->subject, o->firstRun + ir, run_type, it_ini, it_fin, fsmp);

         /* Get pulses in run: */
        int *it_p_ini = &(it_pulse_ini[ir*np_run]); /* Indices of trigger-up frames of this run (rel to start of {val}). */
        int *it_p_fin = &(it_pulse_fin[ir*np_run]); /* Indices of trigger-down frames of this run (rel to start of {val}). */
        
        nes_print_run_pulses(stderr, np_run, it_p_ini, it_p_fin, it_ini, it_fin, fsmp);
        
        assert(it_p_ini[0] >= it_ini);
        assert(it_p_fin[np_run-1] <= it_fin);
        
        /* Decide the marked frame ranges in run: */
        int it_m_ini[nm]; /* Index of first marked frame (rel to start of {val}) for each marker. */
        int it_m_fin[nm]; /* Index of last marked frame (rel to start of {val}) for each marker. */
        for (im = 0, ic = nc_in;  im < nm; im++, ic++)
          { nes_marker_t *mki = &(o->markers.e[im]);
            int kp = mki->ip; /* Index of referece pulse in run. */
            assert((kp >= 0) && (kp < np_run));
            /* Find reference edge {it_m_ref} of reference pulse: */
            int it_m_ref = (mki->edge == 0 ? it_p_ini[kp] : it_p_fin[kp]);
            /* Find frame range {it_m_ini,it_m_fin} of marker: */
            it_m_ini[im] = it_m_ref + mki->shift; 
            it_m_fin[im] = it_m_ini[im] + mki->nt - 1;
            demand((it_m_ini[im] >= it_ini) && (it_m_fin[im] <= it_fin), "marked phase overflows run");
          }
        
        /* Copy the data of this run to {vru[0..nt_run-1][0..nc_out-1]}, set the phase marker channels: */
        double **vru = notnull(malloc(nt_run*sizeof(double*)), "no mem");
        int kt; /* Frame index in {vru}. */
        for (kt = 0; kt < nt_run; kt++)
          { int it = it_ini + kt; /* Frame index in input file. */
            double *valt = val[it];
            double *vrut = notnull(malloc(nc_out*sizeof(double)), "no mem");
            vru[kt] = vrut;
            /* Copy original channels: */
            for (ic = 0; ic < nc_in; ic++) { vrut[ic] = valt[ic]; }
            /* Set the phase marker channels: */
            for (im = 0, ic = nc_in;  im < nm; im++, ic++)
              { nes_marker_t *mki = &(o->markers.e[im]);
                vrut[ic] = ((it >= it_m_ini[im]) && (it <= it_m_fin[im]) ? mki->val : 0.0);
              }
          }
        
        /* Set the header fields specific to this run: */ 
        h->orig->run = o->firstRun + ir;
        h->orig->it_ini = it_ini;
        h->orig->it_fin = it_fin;
        
        h->type = run_type;
        h->nt = nt_run;

        /* Write the selected frames: */ 
        nes_write_eeg_run(o->outDir, h, nt_run, nc_out, vru, np_run, it_p_ini, it_p_fin);
        
        /* Release the storage: */ 
        for (kt = 0; kt < nt_run; kt++) { free(vru[kt]); } 
        free(vru);
      }
    
    /* Don't bother freeing the allocated storage. */
    return 0;
  }

void nes_print_run_pulses(FILE *wr, int np, int it_p_ini[], int it_p_fin[], int it_ini, int it_fin, double fsmp)
  {
    if (np > 0)
      { /* Print the trigger pulses in this run: */
        fprintf(wr, "\n");
        int k;
        for (k = 0; k < np; k++) 
          { int nt_pulse = it_p_fin[k] - it_p_ini[k] + 1;
            fprintf(wr, "   trigger pulse %2d:  %8d samples ( %10.4f s )", k, nt_pulse, nt_pulse/fsmp);
            fprintf(wr, "  %8d .. %8d", it_p_ini[k] - it_ini, it_p_fin[k] - it_ini);
            fprintf(wr, " ( %10.4f - %10.4f s )\n", (it_p_ini[k] - 0.5 - it_ini)/fsmp, (it_p_fin[k] + 0.5 - it_ini)/fsmp);
          }
        fprintf(wr, "\n");
      }
  }

void nes_print_run_info(FILE *wr, char *class, int subject, int ir, char *type, int it_ini, int it_fin, double fsmp)
  { 
    int nt = it_fin - it_ini + 1; /* Number of data frames in run. */
    double ts_ini_g = (it_ini - 0.5)/fsmp; /* Start time of run (seconds) from start of orig file. */
    double ts_fin_g = (it_fin + 0.5)/fsmp; /* End time of run (seconds) from start of orig file. */

    fprintf(wr, "subject %3d run %3d type %-8s (%s)", subject, ir, type, class);
    fprintf(wr, " - %8d samples ( %10.4f s )", nt, nt/fsmp);
    fprintf(wr, " %8d .. %8d", it_ini, it_fin);
    fprintf(wr, " ( %10.4f _ %10.4f s )\n", ts_ini_g, ts_fin_g);
  }
    
void nes_write_eeg_run(char *outDir, neuromat_eeg_header_t *h, int nt, int nc, double **vru, int np, int it_p_ini[], int it_p_fin[])
  {
    int it_ini = h->orig->it_ini;
    int it_fin = h->orig->it_fin; 
    demand((0 <= it_ini) && (it_ini <= it_fin) && (it_fin < h->orig->nt), "invalid frame index range");
    demand(h->nt == nt, "inconsistent {nt} in header");
    demand(h->nc == nc, "inconsistent {nc} in header");
    demand((0 <= h->ne) && (h->ne <= nc), "invalid {ne} in header");
    
    /* Currently we do not downsample, so: */
    demand((! isnan(h->fsmp)) && (h->fsmp == h->orig->fsmp), "inconsistent sampling freqs");
    
    char *fname = NULL;
    asprintf(&fname, "%s/s%03d_r%03d.txt", outDir, h->orig->subject, h->orig->run);
    FILE *wr = open_write(fname, TRUE);
    
    neuromat_eeg_header_write(wr, h);
    
    neuromat_eeg_data_write(wr, nt, nc, vru, 0, nt - 1, 1);
    fclose(wr);
    free(fname);
  }
  
#define nes_MAX_SUBJECT 999
#define nes_MAX_RUN 999
#define nes_MAX_PULSES_PER_RUN 100
#define nes_MAX_FRAMES_PER_RUN 72000
#define nes_MIN_FSMP 0.01
#define nes_MAX_FSMP 1000.0
#define nes_MIN_MARKER_VAL 0.01
#define nes_MAX_MARKER_VAL 1000.0

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
    
    argparser_get_keyword(pp, "-subject");
    o->subject = (int)argparser_get_next_int(pp, 1, nes_MAX_SUBJECT);
    
    if (argparser_keyword_present(pp, "-firstRun"))
      { o->firstRun = (int)argparser_get_next_int(pp, 1, nes_MAX_RUN); }
    else
      { o->firstRun = 1; }
      
    /* Parse channel name list: */
    argparser_get_keyword(pp, "-channels");
    o->chname = nes_parse_next_token_args(pp);
    int nc = o->chname.ne;
    if (nc < 1) { argparser_error(pp, "no channels?"); }
    
    /* Parse trigger channel name option and look it up in channel name list: */
    argparser_get_keyword(pp, "-trigger");
    char *trigger_name = argparser_get_next_non_keyword(pp);
    int ict = 0;
    while ((ict < nc) && (strcmp(trigger_name, o->chname.e[ict]) != 0)) { ict++; }
    if (ict >= nc) { argparser_error(pp, "invalid trigger channel name"); }
    o->trigger_ic = ict;
  
    if (argparser_keyword_present(pp, "-nElectrodes"))
      { o->nElectrodes = (int)argparser_get_next_int(pp, 1, o->trigger_ic); }
    else
      { o->nElectrodes = o->trigger_ic; }
      
    if (argparser_keyword_present(pp, "-runTypes"))
      { o->runTypes = nes_parse_next_token_args(pp); }
    else
      { o->runTypes = string_vec_new(0); }
    
    argparser_get_keyword(pp, "-pulsesPerRun");
    o->pulsesPerRun = (int)argparser_get_next_int(pp, 1, nes_MAX_PULSES_PER_RUN);

    argparser_get_keyword(pp, "-fSampling");
    o->fSampling = argparser_get_next_double(pp, nes_MIN_FSMP, nes_MAX_FSMP);

    if (argparser_keyword_present(pp, "-framesPerRun"))
      { o->framesPerRun = (int)argparser_get_next_int(pp, 1, nes_MAX_FRAMES_PER_RUN); }
    else
      { o->framesPerRun = 0; }
      
    argparser_get_keyword(pp, "-source");
    o->source = argparser_get_next_non_keyword(pp);
      
    argparser_get_keyword(pp, "-outDir");
    o->outDir = argparser_get_next_non_keyword(pp);
    
    /* Parse marker list: */
    o->markers = nes_marker_vec_new(10);
    int nm = 0;
    while(argparser_keyword_present(pp, "-marker"))
      { char *name = argparser_get_next_non_keyword(pp);
        int ip = (int)argparser_get_next_int(pp, 0, o->pulsesPerRun-1); 
        int edge = (int)argparser_get_next_int(pp, 0, 1); 
        int shift = (int)argparser_get_next_int(pp, -(nes_MAX_FRAMES_PER_RUN-1), nes_MAX_FRAMES_PER_RUN-1);
        int nt = (int)argparser_get_next_int(pp, 1, nes_MAX_FRAMES_PER_RUN);
        double val = argparser_get_next_double(pp, nes_MIN_MARKER_VAL, nes_MAX_MARKER_VAL);
        nes_marker_vec_expand(&(o->markers),nm);
        o->markers.e[nm] = 
          (nes_marker_t) { .name = name, .ip = ip, .edge = edge, .shift = shift, .nt = nt, .val = val };
        nm++;
      }
    nes_marker_vec_trim(&(o->markers),nm);

    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }
  
string_vec_t nes_parse_next_token_args(argparser_t *pp)
  { int nc = 0;
    string_vec_t toks = string_vec_new(100);
    while (argparser_next_is_non_keyword(pp))
      { string_vec_expand(&toks,nc);
        toks.e[nc] = argparser_get_next_non_keyword(pp);
        nc++;
      }
    if (nc == 0) { argparser_error(pp, "empty token list"); }
    string_vec_trim(&toks,nc);
    return toks;
  }

vec_typeimpl(nes_marker_vec_t,nes_marker_vec,nes_marker_t);
