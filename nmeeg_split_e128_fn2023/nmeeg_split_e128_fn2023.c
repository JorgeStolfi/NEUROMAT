#define PROG_NAME "nmeeg_split_e128_fn2023"
#define PROG_DESC "Splits a multiple-run 128-electrode EEG data file into individual runs."
#define PROG_VERS "2023-10-27"

#define nmeeg_split_e128_fn2023_C_COPYRIGHT \
  "Copyright © 2023 by the State University of Campinas (UNICAMP)"

/* Last edited on 2023-12-05 23:00:41 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -inDir {IN_DIR} \\\n" \
  "    -subject {SUBJ} \\\n" \
  "    -framesPerSeg {NFR_PRE} {NFR_POS} \\\n" \
  "    -outDir {OUT_DIR} \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program reads from file {IN_FILE} a 128-electrode EEG data file" \
  " from the \"Renewal Points\" experiment by Najman et al [FN2023], and" \
  " extracts from it the segments that correspond to individual stimuli" \
  " as separate files.  It also provides a new `stimulus' marker channel.\n" \
  "\n" \
  "INPUT EEG FILE FORMAT\n" \
  "  The main input file is \"{IN_DIR}/v{VV}/eeg.txt\", where {VV} is the" \
  " two-digit {SUBJ} idenfifier. It should be a full-session EEG recording" \
  " of the 128-electrode equipment, converted to \".txt\" format.  It should" \
  " contains one data frame per lines, each with one sample" \
  " from {NC=NE=128} data channels (electrode potentials in microvolts).\n" \
  "\n" \
  "  The program splits the file into segments, one for each stimulus" \
  " presented. See the  INPUT STIMULUS DATA FILE section below.  Each" \
  " segment will comprise {NFR_PRE} frames before" \
  " the stimulus start time and continues for an additional {NFR_POS} frames" \
  " from that time onwards.\n" \
  "\n" \
  "  The input file should have a header compatible" \
  " with {neuromat_eeg_header_read} in file {neuromat_eeg_header.h}, that " \
  " specifies the sampling frequency and the names of the channels.\n" \
  "\n" \
  "INPUT STIMULUS DATA FILE\n" \
  "  The program also reads a separate file \"{IN_DIR}/v{VV}/sti.txt\" with" \
  " the starting time and type for the stimulus in each segment.  This file" \
  " should have one line per segment, with the segment" \
  " number {NNNN} (sequential, starting from 0), the stimulus start" \
  " time {ST_TIME} (in seconds, from the start of the file), and the" \
  " stimulus type {ST_TYPE} (0, 1, or 2).\n" \
  "\n" \
  "OUTPUT FILES\n" \
  "  The output files are written with" \
  " names \"{OUT_DIR}/v{VV}/n{NNNN}.txt\" where {VV}" \
  " is the two-digit subject ID number (obtained from the input file's header), {NNNN} is" \
  " a four-digit sequential segment number within the file. The" \
  " directory \"{OUT_DIR}/v{VV}\" must exist. \n" \
  "\n" \
  "  Each output file contains the specified data frames for a single" \
  " segment corresponding to a single stimulus.\n" \
  "\n" \
  "  Each output frame contains {NE+1} channels, comprising the same {NE} electrode" \
  " channels as in the input (including the \"CZ\" electrode) followed by" \
  " an additional marker channel \"ST\".  The value of the latter is the stimulus type -- 0, 1, or 2.\n" \
  "\n" \
  "  Each output file begins with header records that specify some relevant" \
  " parameters such as number of frames and channels, sampling frequency, etc..\n" \
  "\n" \
  "OPTIONS\n" \
  "  -inDir {IN_DIR}\n" \
  "    This mandatory argument specifies the directory containing the EEG and segment data file.\n" \
  "\n" \
  "  -subject {SUBJ}\n" \
  "    This mandatory argument specifies the ID number of the subject" \
  " that yielded the input EEG file.  Should be a positive number.\n" \
  "\n" \
  "  -framesPerSeg {NFR_PRE} {NFR_POS}\n" \
  "    This mandatory argument specifies the number of frames to extract" \
  " for each segment; specifically, {NFR_PRE} frames before the stimulus start time" \
  " plus {NFR_POS} frames starting at that time.\n" \
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
  "  Created 2023-10-27 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2023-10-27 Created based on {nmeeg_split_e128_gh2013.c}. (J.Stolfi).\n" \
  "  2023-11-01 Assuming the input file has a NEUROMAT standard header.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " nmeeg_split_e128_fn2023_C_COPYRIGHT ".\n" \
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
#include <fget.h>
#include <affirm.h>
#include <jsfile.h>
#include <jsmath.h>
#include <jsstring.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_io.h>
#include <neuromat_eeg_header.h>
#include <neuromat_eeg_source.h>
#include <neuromat_eeg_frame_buffer.h>

typedef struct nesf_options_t
  { 
    char *inDir;              /* INput file name. */
    int32_t subject;           /* Subject ID number. */
    int32_t framesPerSeg_pre;  /* Desired number of data frames before start-stimulus frame. */
    int32_t framesPerSeg_pos;  /* Desired number of data frames from start-stimulus frame. */
    char *outDir;              /* Directory for all output files. */
  } nesf_options_t;
  /* Arguments from command line. */
  
#define nesf_MAX_SUBJECT 99
  /* Max subject ID number. */

#define nesf_MAX_SEGMENT 9999
  /* Max segment ID number within an experimental block. */

#define nesf_MAX_FRAMES_PER_SEGMENT 72000
  /* Max frames in each output segment that the user may specify
    in the command line. */

char *nesf_input_file_name(char *dir, int32_t subject, char *name);
  /* Returns the string "{dir}/v{VV}/name.txt", where {VV} is the 
    two-digit {subject} id number. */

neuromat_eeg_header_t *nesf_read_header(FILE *rd, int *nl_P);
  /* Reads from {rd} a neuroMat EEG dataset header, assuming the ".txt" 
    format, such as created by {neuromat_eeg_header_write}.
    
    If {nl_P} is not null, {*nl_P} is assumed to be the number of
    file lines (including headers and comments) previously read from
    {rd}; usually 0. The count is updated if not null. */

bool_t nesf_read_segment_data
  ( FILE *rd,
    int32_t *idSeg_P,
    double *stTime_P,
    int32_t *stType_P,
    int32_t *nl_P,
    int32_t *ns_P
  );
  /* Tries to read from {rd} the next data line containing the ID number
    of a segment (a non-negative integer), the time of the start of the
    stimulus (in seconds, from the first input frame time), and the
    stimulus type (0, 1, or 2). If there is such a line, stored those
    values in {*idSeg_P,*stTime_P,*stType_P} and return {TRUE}. 
    
    If there is no next data line, returns {FALSE}.
    
    Ignores '#'-comments and skips lines with no data. Increments {*nl_P)
    for each line read, and {*ns_P} if it finds a data line. */

double **nesf_alloc_segment(int32_t nt, int32_t nc);
  /* Allocate an EEG dataset {vru} with {nt} frames and {nc} channels,
    specifically {vru[0..nt-1][0..nc-1]}. */
        
neuromat_eeg_header_t *nesf_make_segment_header(neuromat_eeg_header_t *hin);
  /* Creates an EEG dataset header {hot} for a generic output segment, compatible with
    {neuromat_eeg_header_write}.  Most data fields are copied from the input header
    {hin}.  The fields {hot->{type,nt,run}}, {hot->orig->{it_ini,it->fin}} are
    left undefined. */

void nesf_extract_segment
  ( neuromat_eeg_frame_buffer_t *buf,
    neuromat_eeg_frame_buffer_read_proc_t *read_frame, 
    int32_t it_ini,
    int32_t it_fin, 
    int32_t nt,
    int32_t nc,
    double **vru
  );
  /* Extracts a single experiment segment from the EEG dataset buffer {buf},
    comprising frames {it_ini..it_fin}, and stores them in
    {vru[0..nt-1][0..nc-1]}.
    
    The procedure assumes that the values for each frame with index {it}
    in that interval are in {buf.val[itb][0..buf.nc-1]} where
    {itb=(it%buf.size)}. It also requires {nt==it_fin-it_ini+1} and
    {nc==buf.nc}. The array {vru} must have been allocated by the
    client.
    
    Frame {it_ini} must not have been flushed from buffer {buf} yet. Reads
    more frames from {rd} if needed with the client-given procedure
    {read_frame}, which should return 1 if it succeeds, 0 if there
    are no more frames.  
    
    Fails if {it_ini} is negative or {it_fin} is past the last frame. */

void nesf_add_marker_channels(int32_t nt, int32_t nc, double **vru, int32_t nm, int32_t stType);
  /* Adds {nm} new /marker/ channels to the EEG segment data {vru}.
  
    Assumes that {vru} has been allocated with {nt} rows and {nc+nm}
    columns. Add the {nm} new channels as columns {nc+k} where {k}
    ranges in {0..nm-1}. Currently {nm} must be 1, and the new maker
    channel ("ST", the stimulus type) has constant value {stType} on all
    frames of the segment. */

void nesf_adjust_electrode_averages_in_segment
  ( int32_t nt,
    int32_t nc,
    double **vru,
    int32_t nt_pre,
    int32_t ne,
    bool_t verbose
  );
  /* Shifts all electrode channels {0..ne-1} in the frames
    {vru[0..nt-1]} of a single segment so that the mean value of
    the first {nt_ini} samples is zero. Channels {nc..ne-1}, if any,
    are not modified.  Ignored apparent outliers. */

double nesf_adjust_electrode_average
  ( int32_t nt,
    int32_t nc,
    double **vru,
    int32_t ie,
    double wt[], /* (IN) Frame weights. */
    double pr[]  /* (WORK) Validity prob. */
  );
  /* Shifts the samples {vru[0..nt-1][ie]} of electrode {ie} so that the value of
    each electrode is zero.  In the mean computation, each frame {vru[kt]} 
    gets weight {wt[kt]}, and is further un-wheighted if its sample
    appears to be an outlier.  The scratch array {pr} must have {nt} elements. */
          
double nesf_bayes(double v, double pri_gud, double avg_gud, double dev_gud, double avg_bad, double dev_bad);
  /* Uses Bayes's formula to compute the probability of {v} being a good sample,
    assuming good and bad samples have normal distributions with averages
    {avg_gud,avg_bad} and deviations {dev_gud,dev_bad}, respectively, 
    and that, a priori, the sample is good with probability {pri_gud}. */

void nesf_compute_stats
  ( int32_t nt, 
    int32_t nc, 
    double **vru, 
    int32_t ie, 
    double wt[], 
    double pr[], 
    bool_t good, 
    double *avg_P, 
    double *dev_P, 
    double *pri_P
  );
  /* Estimates the average {*avg_P} and deviation {*dev_P} of samples {vru[0..nt-1][ie]}.
    If {good} is TRUE, estimates the parameters of the `good' samples (inliers),
    if {good} is FALSE estimates those of the `bad' ones (outliers).
    Also returns in {*pri_P} (if not null) the estimated `a priori' probability of a 
    sample being good or bad, respectively.

    Assumes that {wt[kt]} is the relevance weight of sample {vru[kt][ie]},
    and {pr[kt]} is the current probability of that sample being a good sample. 
    If {pr} is null, assumes {pr[kt]==0.5} for all {kt}. */

void nesf_compute_frame_weights_in_segment
  ( int32_t nt,
    int32_t nc,
    double **vru,
    int32_t nt_pre,
    double *wt
  );
  /* Computes sample averaging weights for each frame in the segment.
    Assumes that the frames are {vru[0..nt-1][0..nc-1]}. Currently, the weight is
    1 for frames {0..nt_pre-1} and 0 for the rest. */

void nesf_output_segment
  ( int32_t nt,
    int32_t nc,
    double **vru,
    char *outDir,
    neuromat_eeg_header_t *hot
  );
  /* Writes to disk one segment extracted from the EEG dataset, comprising
    frames frames {vru[0..nt-1][0..nc-1]}.
    
    The output file will be called "{outDir}/v{VV}/n{NNNN}.txt", where
    {VV} is the two-digit value of {hot.subject}, and {NNNN} is
    the four-digit value of {hot.run}, both zero-padded. The
    directory "{outDir}/v{VV}/" must exist.
    
    The procedure sets {hot.nt=nt} and {hot.nc=nc}, and assumes that the client has
    defined all other fields of {hot} properly, including {hot{.ne,.type,.subject,.run,.chname[0..nc-1]}},
    {hot.orig.{it_ini,.it_fin}}. Requires {hot.subject == hot.orig.subject}
    and {hot.orig.run == INT32_MIN}.
    
    The procedure then uses {neuromat_eeg_header_write} to write the
    non-null information from {hot}, then calls {neuromat_eeg_data_write}
    to write the data. */

void nesf_print_segment_info(FILE *wr, char *class, neuromat_eeg_header_t *hot);
  /* Prints to {wr} a one-line summary for one experiment segment. Assumes that the segment has the
    given {hot.type} and spans from frame {hot.orig.it_ini} to frame {hot.orig.it_fin} of the
    original data file. The string {class} is printed in the summary.
    The sampling frequency {hot.fsmp} is used to convert frame indices into
    time coordinates. */

void nesf_print_segment_phases
  ( FILE *wr, 
    int32_t it_ini, 
    int32_t it_fin, 
    int32_t nt_pre,
    char *name[],
    double fsmp
  );
  /* Prints to {wr} the start and stop of the two phases (preamble and response)
    of a segment that spans frames with indices {it_ini..it_fin} of the original data file.
    Assumes that the premble comprises the first {nt_pre} frames.  The sampling frequency
    {fsmp} is used to convert frame indices into time coordinates. */
      
nesf_options_t *nesf_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments. */
  
int32_t main(int32_t argc, char **argv);
  /* Main prog. */

int32_t main(int32_t argc, char **argv)
  {
    nesf_options_t *o = nesf_parse_options(argc, argv);
    bool_t debug = TRUE;
    bool_t debugFrame = FALSE; /* Print debug for each frame read. */

    bool_t verbose = TRUE;
    
    /* Print some options: */
    fprintf(stderr, "extracting %d data frames before stimulus onset\n", o->framesPerSeg_pre);
    fprintf(stderr, "extracting %d data frames on and after stimulus onset\n", o->framesPerSeg_pos);

    fprintf(stderr, "output files are named \"%s/v%02d/n{NNNN}.txt\"\n", o->outDir, o->subject);
    
    /* Input EEG file: */
    char *fname_eeg = nesf_input_file_name(o->inDir, o->subject, "eeg");
    FILE *rd_eeg = open_read(fname_eeg, verbose);
    int32_t nl_eeg = 0; /* Number of lines read so far from {rd_eeg}. */
    int32_t nt_eeg = 0; /* Number of data frames read so far from {rd_eeg}. */

    /* Read the file header: */
    neuromat_eeg_header_t *hin = nesf_read_header(rd_eeg, &nl_eeg);

    int32_t nt = hin->nt; /* Count of data frames expected in input file. */
    int32_t nc_in = hin->nc; /* Total count of data channels (incl. vref, triggers). */
    int32_t ne = hin->ne; /* Count of electrode channels (incl. vref). */
    
    int32_t nm = 1; /* Number of extra marker channels to add to output files. */
    int32_t nc_out = ne + nm; /* Number of channels in output files: */

    /* Some data for the benefit of scripts that parse the {stderr} output: */
    fprintf(stderr, "subject = %02d # Subject ID.\n", hin->subject);
    fprintf(stderr, "nc_in = %d     # Count of input channels.\n", nc_in);
    fprintf(stderr, "ne = %d        # Count of electrode channels.\n", ne);
    fprintf(stderr, "fsmp = %f      # Sampling frequency.\n", hin->fsmp);
    fprintf(stderr, "nm = %d        # Count of added maker channels.\n", nm);
    fprintf(stderr, "nc_out = %d    # Count of output channels.\n", nc_out);

    /* Allocate the extracted segment frames: */
    int32_t nt_segment = o->framesPerSeg_pre + o->framesPerSeg_pos; /* Frames per segment. */
    double **vru = nesf_alloc_segment(nt_segment, nc_in);
    /* Allocate the EEG data frame buffer: */
    int32_t nt_buf = 2*nt_segment; /* Frames to keep in buffer; should be enough. */
    neuromat_eeg_frame_buffer_t *buf = neuromat_eeg_frame_buffer_new(nt_buf, nc_in);

    /* Header for output files: */
    neuromat_eeg_header_t *hot = nesf_make_segment_header(hin);
    assert(hot->nc == nc_out);
    assert(hot->ne == ne);
    assert(hot->orig->run == INT32_MIN);
    
    /* Stimulus data file: */
    char *fname_sti = nesf_input_file_name(o->inDir, o->subject, "sti");
    FILE *rd_sti = open_read(fname_sti, verbose);
    int32_t nl_sti = 0; /* Number of lines read so far from {rd_sti}. */
    int32_t ns_sti = 0; /* Number of stimulus data lines read so far from {rd_sti}. */

    auto int32_t read_frame(int32_t nc_fr, double val_fr[]);
      /* Reads a frame from {rd_eeg} (in the plain text format)
        into {val_fr[0..nc_fr-1]}.  Expects {nc_fr==nc_in}. 
        Increments {nl_eeg} and {nt_eeg} as appropriate. */

    /* Main loop: Read segments from {rd_sti,rd_eeg}, write them out. */
    int32_t ns = 0;  /* Number of segments processed. */
    int32_t it_next; /* Index of next un-extracted frame in file. */
    while (TRUE)
      { 
        /* Try to read segment data from segment file: */
        int32_t idSeg, stType;
        double stTime;
        bool_t ok = nesf_read_segment_data(rd_sti, &idSeg, &stTime, &stType, &nl_sti, &ns_sti);
        if (! ok) { break; }

        fprintf(stderr, "---------------------------------------------------\n");
        if (debug) { fprintf(stderr, "trying to read segment v%02d time = %.6f type = %d\n", idSeg, stTime, stType); }
        demand((idSeg >= 0) & (idSeg <= nesf_MAX_SEGMENT), "invalid segment number");
        
        /* Read segment from {rd_eeg}: */
        int32_t it_sti = (int32_t)floor(stTime*hin->fsmp + 0.5); /* Start-of-stimulus frame. */
        int32_t it_ini = it_sti - o->framesPerSeg_pre;
        int32_t it_fin = it_sti + o->framesPerSeg_pos - 1;
        if (it_ini > it_next)
          { fprintf(stderr, "skipping %d frames {%d .. %d}\n", it_ini - it_next, it_next, it_ini-1); }
        fprintf(stderr, "frames %d..%d, stimuls start at %d\n", it_ini, it_fin, it_sti);
        demand(it_ini >= it_next, "segment starts at negative frame index or overlaps previous one");
        
        /* Extract the segment frames to {vru[0..nt-1][o..nc_in-1]}: */
        nesf_extract_segment(buf, read_frame, it_ini, it_fin, nt_segment, nc_in, vru);

        /* Redefine the marker channels: */
        nesf_add_marker_channels(nt_segment, nc_in, vru, nm, stType);
        
        /* Shift electrode averages to zero: */
        bool_t show_avg = verbose;
        nesf_adjust_electrode_averages_in_segment(nt_segment, nc_out, vru, o->framesPerSeg_pre, ne, show_avg);

        /* Write segment to disk: */
        hot->nt = nt_segment;
        hot->subject = hin->subject;
        hot->run = idSeg;
        hot->type = fmt_int(stType, 1);
        
        hot->orig->it_ini = it_ini;
        hot->orig->it_fin = it_fin;
        
        nesf_output_segment(nt_segment, nc_out, vru, o->outDir, hot);
        
        if (debug) { fprintf(stderr, "done segment v%02d\n", idSeg); }

        /* Count one more segment, prepare for next one: */
        ns++;

        /* Assume that pulses do not overlap: */
        it_next = it_fin + 1;
        fprintf(stderr, "---------------------------------------------------\n");
      }
      
    /* Read and ignore rest of EEG file: */
    if (it_next < nt)
      { fprintf(stderr, "ignoring %d final frames {%d .. %d}\n", nt-it_next, it_next, nt-1); }
      
    fclose(rd_eeg);
    fclose(rd_sti);
        
    /* Final reports: */
    fprintf(stderr, "read %d EEG data frames in %d file lines\n", nt_eeg, nl_eeg);
    demand(nt_eeg == nt, "frame count in file does not match count in header");
    fprintf(stderr, "read %d stimulus data lines in %d file lines\n", ns_sti, nl_sti);
    demand(ns > 0, "no segments detected in input file");
    demand(ns_sti == ns, "stimulus count in file does not match output segment count");
   
    /* Release some storage: */ 
    fprintf(stderr, "freeing {vru}\n");
    for (int32_t kt = 0; kt < nt_segment; kt++) { free(vru[kt]); } 
    free(vru);
    neuromat_eeg_frame_buffer_free(buf);
    
    return 0;
    
    /* Internal implementations: */
    int32_t read_frame(int32_t nc_fr, double val_fr[])
      { assert(nc_fr == buf->nc);
        if (debugFrame) { fprintf(stderr, "  trying to read frame %d at line %d\n", nt_eeg+1, nl_eeg+1); }
        int32_t nrd = neuromat_eeg_frame_read(rd_eeg, buf->nc, val_fr, &nl_eeg, &nt_eeg);
        assert((nrd == 0) || (nrd == 1));
        return nrd;
      }
  }

char *nesf_input_file_name(char *dir, int32_t subject, char *name)
  { char *fname = NULL;
    asprintf(&fname, "%s/v%02d/%s.txt", dir, subject, name);
    return fname;
  }

bool_t nesf_read_segment_data(FILE *rd, int32_t *idSeg_P, double *stTime_P, int32_t *stType_P, int32_t *nl_P, int32_t *ns_P)
  { while (TRUE)
      { if (fget_test_comment_or_eol(rd, '#', NULL)) { (*nl_P)++; continue; }
        if (fget_test_eof(rd)) { return FALSE; }
        break;
      }
    (*nl_P)++;
    int32_t idSeg = fget_int32(rd);
    double stTime = fget_double(rd);
    int32_t stType = fget_int32(rd);
    fget_comment_or_eol(rd, '#', NULL);
    demand((idSeg >= 0) && (idSeg <= nesf_MAX_SEGMENT), "invalid segment number");
    demand(idSeg == (*ns_P), "segment numbers out of sequence");
    demand(stTime >= 0.0, "invalid stimulus start time");
    demand((stType >= 0) && (stType <= 2), "invalid stimulus type");
    (*idSeg_P) = idSeg;
    (*stTime_P) = stTime;
    (*stType_P) = stType;
    (*ns_P)++;
    return TRUE;
  }

neuromat_eeg_header_t *nesf_make_segment_header(neuromat_eeg_header_t *hin)
  {
    neuromat_eeg_header_t *hot = neuromat_eeg_header_new();
    neuromat_eeg_header_merge(hot, hin);
    neuromat_eeg_header_append_marker_channel(hot, "ST");
    assert(hin->orig->subject == hin->subject);
    assert(hin->run == INT32_MIN); /* All runs. */
    
    /* Original file data (for merge purpose): */
    hot->orig = neuromat_eeg_source_new();
    neuromat_eeg_header_merge_orig(hot->orig, hin->orig);
    assert(hot->orig->run == INT32_MIN); /* All runs. */
    assert(hot->orig->subject == hin->subject);
    
    /* These fields are redefined later for each segment: */
    hot->nt = INT32_MIN;
    hot->type = NULL;
    hot->orig->it_ini = INT32_MIN; 
    hot->orig->it_fin = INT32_MIN;

    return hot;
  }

neuromat_eeg_header_t *nesf_read_header(FILE *rd, int *nl_P)   
  {
    neuromat_eeg_header_t *h = neuromat_eeg_header_read(stdin, 128, 500, nl_P);

    fprintf(stderr, "original file name = \"%s\"\n", h->orig->file);
    fprintf(stderr, "subject ID = %03d\n", h->subject);
    fprintf(stderr, "sampling frequency = %.10g\n", h->fsmp);

    fprintf(stderr, "expects %d data frames\n", h->nt);
    fprintf(stderr, "each input frame has %d channels =", h->nc);
    for (int ic = 0; ic < h->nc; ic++) 
      { if ((ic % 10) == 0) { fprintf(stderr, "\n "); }
        fprintf(stderr, " %s", h->chname[ic]);
      }
    fprintf(stderr, "\n");
    /* Include the "CZ" channel in the electrode count, if present: */
    if ((h->ne < h->nc) && (strcmp(h->chname[h->ne], "CZ") == 0))
      { fprintf(stderr, "including \"CZ\" as an electrode\n");
        h->ne++;
      }
    fprintf(stderr, "channels 0..%d (%d channels) are electrodes\n", h->ne-1, h->ne);
    return h;
  }


double **nesf_alloc_segment(int32_t nt, int32_t nc)
  {
    double **vru = talloc(nt, double*);
    for (int32_t it = 0; it < nt; it++) { vru[it] = talloc(nc, double); }
    return vru;
  }

void nesf_extract_segment
  ( neuromat_eeg_frame_buffer_t *buf,
    neuromat_eeg_frame_buffer_read_proc_t *read_frame,
    int32_t it_ini,
    int32_t it_fin,
    int32_t nt,
    int32_t nc,
    double **vru
  )   
  {
    assert((nc > 0) && (nc == buf->nc));
    assert((nt > 0) && (nt == it_fin - it_ini + 1));
    
    /* Make sure frames {it_ini..it_fin} are in buffer: */
    demand(nt <= buf->size, "not enough slots in buffer");
    bool_t mirror = FALSE; /* Should generate missing frames by mirroring? */


    for (int32_t it = it_ini; it <= it_fin; it++)
      { int32_t kt = it - it_ini; /* Index of frame {it} in {vru}. */
        assert(kt < nt);
        /* Read frame {it}.  Let the routine handle mirroring: */
        int32_t itb = neuromat_eeg_frame_buffer_get_frame(buf, it, read_frame, mirror);
        demand(itb >= 0, "frame not found in file");
        assert(itb < buf->size);
        int32_t ic;
        for (ic = 0; ic < nc; ic++) { vru[kt][ic] = buf->val[itb][ic]; }
      }
  }

;
  /* Adds {nm} new /marker/ channels to the EEG segment data {vru}.
  
    Assumes that {vru} has been allocated with {nt} rows and {nc+nm}
    columns. Add the {nm} new channels as columns {nc+k} where {k}
    ranges in {0..nm-1}. Currently {nm} must be 1, and the new maker
    channel ("ST", the stimulus type) has constant value {stType} on all
    frames of the segment. */

void nesf_add_marker_channels(int32_t nt, int32_t nc, double **vru, int32_t nm, int32_t stType)
  { demand(nm == 1, "Wrong number of marker channels");
    int ncols = nc + nm;
    for (int32_t kt = 0; kt < nt; kt++)
      { vru[kt][ncols-1] = (double)stType; }
  }

void nesf_adjust_electrode_averages_in_segment
  ( int32_t nt,
    int32_t nc,
    double **vru,
    int32_t nt_pre,
    int32_t ne,
    bool_t verbose
  )     
  {
    double *wt = talloc(nt, double);  /* Importance weights. */
    double *pr = talloc(nt, double);  /* Validity probs. */
    nesf_compute_frame_weights_in_segment(nt, nc, vru, nt_pre, wt);
    for (int32_t ie = 0; ie < ne; ie++)
      { if (verbose) { fprintf(stderr, "    electrode %3d", ie); }
        double avg = nesf_adjust_electrode_average(nt, nc, vru, ie, wt, pr);
        if (verbose) { fprintf(stderr, " avg %+12.4f\n", avg); }
      }
    free(wt);
    free(pr);
  }
    
void nesf_compute_frame_weights_in_segment
  ( int32_t nt,
    int32_t nc,
    double **vru,
    int32_t nt_pre,
    double *wt
  )
  { 
    /* Assign weight 1 to first {nt_pre} frames, 0 for the rest: */
    for (int32_t kt = 0; kt < nt; kt++) { wt[kt] = (kt < nt_pre ? 1.0 : 0.0);  }
    /* !!! Perhaps de-emphasize near edges? */
  }
   
double nesf_adjust_electrode_average
  ( int32_t nt,
    int32_t nc,
    double **vru,
    int32_t ie,
    double wt[],
    double pr[]
  )     
  {
    bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "\n  --- nesf_adjust_electrode_average ---\n"); }
    demand((ie >= 0) && (ie < nc), "invalid channel index");
    
    /* Paramters for the overall distribution: */
    double avg_all = NAN;   /* Assumed average of true samples. */
    double dev_all = NAN;   /* Assumed deviation of true samples. */
    nesf_compute_stats(nt, nc, vru, ie, wt, NULL, TRUE,  &avg_all, &dev_all, NULL);
    assert(! isnan(avg_all));
    assert(! isnan(dev_all));
    
    double avg_use;  /* Average value to subtract from all samples. */
    
    if (dev_all < 1.0e-6)
      { /* Signal is dead. */
        avg_use = avg_all; 
      }
    else
      { /* Estimate the average excluding outliers: */
    
        /* Parameters of good samples (updated at each iteration): */
        double avg_gud = avg_all;      /* Assumed average of good samples. */
        double dev_gud = 0.5*dev_all;  /* Assumed deviation of good samples. */
        double pri_gud = 0.5;          /* A priori probability of sample being good. */

        /* Parameters of outlier distribution (maybe to be adjusted): */
        double avg_bad = avg_all;      /* Assumed average of bad samples. */
        double dev_bad = 8.0*dev_all;  /* Assumed deviation of bad samples. */
        double pri_bad = 1 - pri_gud;  /* A priori probability of sample being bad. */

        bool_t adjust_bad_distr = FALSE; /* TRUE to adjust the outlier distribution too. */

        int32_t max_iters = 5; /* Max iterations (was 5). */
        for (int32_t iter = 0; iter < max_iters; iter++)
          { if (debug) { fprintf(stderr, "    iteration %d\n", iter); }
            /* Cook the deviations if the distributions are equal: */
            bool_t same_avg = (fabs(avg_gud - avg_bad) < 0.05*dev_all);
            bool_t same_dev = (fabs(dev_gud - dev_bad) < 0.05*dev_all);
            if (same_avg && same_dev) { dev_gud = dev_gud/1.5; dev_bad = dev_bad*1.5; }
            if (debug) { fprintf(stderr, "      gud[0]:  avg = %+9.4f  dev = %9.4f  pri = %9.4f\n", avg_gud, dev_gud, pri_gud); }
            if (debug) { fprintf(stderr, "      bad[0]:  avg = %+9.4f  dev = %9.4f  pri = %9.4f\n", avg_bad, dev_bad, pri_bad); }
            /* Decide the probability {pr[kt]} of each sample {vru[kt][ie]} being an inlier: */
            for (int32_t kt = 0; kt < nt; kt++) 
              { pr[kt] = nesf_bayes(vru[kt][ie], pri_gud, avg_gud, dev_gud, avg_bad, dev_bad); }
            /* Re-estimate the 'good' (inlier) and 'bad' (outlier) population parameters: */
            nesf_compute_stats(nt, nc, vru, ie, wt, pr, TRUE,  &avg_gud, &dev_gud, &pri_gud);
            if (adjust_bad_distr)
              { /* Re-estimate the outlier distribution: */
                nesf_compute_stats(nt, nc, vru, ie, wt, pr, FALSE, &avg_bad, &dev_bad, &pri_bad);
              }
            else
              { /* Keep the outlier distribution, just recopute the a priori prob: */
                pri_bad = 1 - pri_gud;
              }
            if (debug) { fprintf(stderr, "      gud[1]:  avg = %+9.4f  dev = %9.4f  pri = %9.4f\n", avg_gud, dev_gud, pri_gud); }
            if (debug) { fprintf(stderr, "      bad[1]:  avg = %+9.4f  dev = %9.4f  pri = %9.4f\n", avg_bad, dev_bad, pri_bad); }
            assert(fabs(pri_bad + pri_gud - 1.0) < 0.0001);
          }
        avg_use = avg_gud;
      }
      
    /* Subtract average from signal: */
    for (int32_t kt = 0; kt < nt; kt++)
      { double *vt = vru[kt];
        vt[ie] -= avg_use;
      }
    if (debug) { fprintf(stderr, "\n"); }
    return avg_use;
  }
        
void nesf_compute_stats
  ( int32_t nt, 
    int32_t nc, 
    double **vru, 
    int32_t ie, 
    double wt[], 
    double pr[], 
    bool_t good, 
    double *avg_P, 
    double *dev_P, 
    double *pri_P
  )
  {
    /* Compute the average and overall probability: */
    double sum_wpv = 0;
    double sum_wp = 0;
    double sum_w = 1.0e-200;
    int32_t kt;
    for (kt = 0; kt < nt; kt++) 
      { double vt = vru[kt][ie];
        double wk = wt[kt];
        double pk = (pr == NULL ? 0.5 : (good ? pr[kt] : 1.0 - pr[kt]));
        sum_wpv += wk*pk*vt;
        sum_wp += wk*pk;
        sum_w += wk;
      }
    double avg = sum_wpv/sum_wp;
    double pri = sum_wp/sum_w;
    /* Compute the standard deviation: */
    double sum_wpd2 = 0;
    sum_wp = 1.0e-200;
    for (kt = 0; kt < nt; kt++) 
      { double dk = vru[kt][ie] - avg;
        double wk = wt[kt];
        double pk = (pr == NULL ? 0.5 : (good ? pr[kt] : 1.0 - pr[kt]));
        sum_wpd2 += wk*pk*dk*dk;
        sum_wp += wk*pk;
      }
    double dev = sqrt(sum_wpd2/sum_wp); /* Should compensate for bias... */
    /* Return: */
    assert(! isnan(avg));
    assert(! isnan(dev));
    assert(! isnan(pri));
    (*avg_P) = avg;
    (*dev_P) = dev;
    if (pri_P != NULL) { (*pri_P) = pri; }
  }

double nesf_bayes(double v, double pri_gud, double avg_gud, double dev_gud, double avg_bad, double dev_bad)
  {
    double eps = 1.0e-12;
    double dk_gud = (v - avg_gud)/dev_gud;
    double dk_bad = (v - avg_bad)/dev_bad;
    if (dk_gud > 6.0)
      { return eps; }
    else if (dk_bad > 6.0)
      { return 1.0 - eps; }
    else
      { /* Fudge a priori probs to avoid zeros: */
        pri_gud = (pri_gud + eps)/(1.0 + 2*eps);
        double pri_bad = 1 - pri_gud;
        double P_gud = eps + exp(-dk_gud*dk_gud/2)/dev_gud/sqrt(2*M_PI); /* Prob of sample, assuming good. */
        double P_bad = eps + exp(-dk_bad*dk_bad/2)/dev_bad/sqrt(2*M_PI); /* Prob of sample, assuming bad. */ 
        double PP_gud = pri_gud*P_gud;
        double PP_bad = pri_bad*P_bad;
        double pos_gud = PP_gud/(PP_gud + PP_bad);
        pos_gud = (pos_gud + eps)/(1.0 + 2*eps);
        return pos_gud;
      }
  }

void nesf_print_segment_info(FILE *wr, char *class, neuromat_eeg_header_t *hot)
  { 
    int32_t it_ini = hot->orig->it_ini;
    int32_t it_fin = hot->orig->it_fin; 
    double fsmp = hot->fsmp;

    int32_t nt = it_fin - it_ini + 1; /* Number of data frames in segment. */
    double ts_ini_g = (it_ini - 0.5)/fsmp; /* Start time of segment (seconds) from start of orig file. */
    double ts_fin_g = (it_fin + 0.5)/fsmp; /* End time of segment (seconds) from start of orig file. */

    fprintf(wr, "subject %3d segment %4d type %-8s (%s)", hot->subject, hot->run, hot->type, class);
    fprintf(wr, " - %8d samples ( %10.4f s )", nt, nt/fsmp);
    fprintf(wr, " %8d .. %8d", it_ini, it_fin);
    fprintf(wr, " ( %10.4f _ %10.4f s )\n", ts_ini_g, ts_fin_g);
  }

void nesf_output_segment
  ( int32_t nt,
    int32_t nc,
    double **vru,
    char *outDir,
    neuromat_eeg_header_t *hot
  )
  {
    /* Set {nc,nt} in header: */
    assert(nt > 0);
    assert(nc > 0);
    hot->nt = nt;
    hot->nc = nc;

    /* Check header data for this segment: */
    assert((!isnan(hot->fsmp)) && (hot->fsmp == hot->orig->fsmp)); /* Data is assumed to be unfiltered. */
    assert((hot->ne > 0) && (hot->ne <= hot->nc-2)); /* There are at least two marker channels. */
    assert((hot->subject > 0) && (hot->subject <= nesf_MAX_SUBJECT));
    assert((hot->subject > 0) && (hot->subject <= nesf_MAX_SUBJECT));
    assert(hot->type != NULL);
    assert(hot->chname != NULL);
    assert(hot->component == NULL);
    assert(hot->kfmax == INT32_MIN); /* The {kfmax} is irrelevant for time-domain data.  */
    assert(isnan(hot->flo0));
    assert(isnan(hot->flo1));
    assert(isnan(hot->fhi1));
    assert(isnan(hot->fhi0));
    assert(hot->finvert == INT32_MIN);
    
    /* Check header data re original file (note: mirroring allows out-of-range indices): */
    assert((hot->orig->file != NULL) && (strlen(hot->orig->file) > 0));
    assert(hot->orig->nt > 0);
    if ((hot->orig->it_ini < 0) || (hot->orig->it_fin >= hot->orig->nt))
      { fprintf(stderr, "!! warning: segment frame range %d..%d", hot->orig->it_ini, hot->orig->it_fin);
        fprintf(stderr, " extends outside raw file range %d..%d\n", 0, hot->orig->nt - 1);
      }
    assert(hot->orig->it_ini <= hot->orig->it_fin);
    assert((!isnan(hot->orig->fsmp)) && (hot->orig->fsmp > 0));
    assert(hot->orig->subject == hot->subject);
    assert((hot->orig->run > 0) && (hot->orig->run <= nesf_MAX_SEGMENT));
    
    nesf_print_segment_info(stderr, "extracted", hot);
    
    char *fname = NULL;
    asprintf(&fname, "%s/v%02d/n%03d.txt", outDir, hot->subject, hot->run);
    FILE *wr = open_write(fname, TRUE);
    
    neuromat_eeg_header_write(wr, hot);
    
    neuromat_eeg_data_write(wr, nt, nc, vru, "%14.8e", 0, nt - 1, 1);
    fclose(wr);
    free(fname);
  }

nesf_options_t *nesf_parse_options(int32_t argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    nesf_options_t *o = notnull(malloc(sizeof(nesf_options_t)), "no mem"); 
    
    /* Parse keyword parameters: */
    
    argparser_get_keyword(pp, "-inDir");
    o->inDir = argparser_get_next_non_keyword(pp);
    
    argparser_get_keyword(pp, "-subject");
    o->subject = (int32_t)argparser_get_next_int(pp, 1, nesf_MAX_SUBJECT);
    
    argparser_get_keyword(pp, "-framesPerSeg");
    o->framesPerSeg_pre = (int32_t)argparser_get_next_int(pp, 1, nesf_MAX_FRAMES_PER_SEGMENT);
    o->framesPerSeg_pos = (int32_t)argparser_get_next_int(pp, 1, nesf_MAX_FRAMES_PER_SEGMENT);
    
    argparser_get_keyword(pp, "-outDir");
    o->outDir = argparser_get_next_non_keyword(pp);
  
    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }
