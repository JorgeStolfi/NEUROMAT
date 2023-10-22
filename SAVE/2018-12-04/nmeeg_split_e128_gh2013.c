#define PROG_NAME "nmeeg_split_e128_gh2013"
#define PROG_DESC "Splits a multiple-run 128-electrode EEG data file into individual runs."
#define PROG_VERS "2013-06-05"

#define nmeeg_split_e128_gh2013_C_COPYRIGHT \
  "Copyright © 2013 by the State University of Campinas (UNICAMP)"

/* Last edited on 2023-10-21 21:55:25 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    [ -firstBlock {FBLOCK} ] \\\n" \
  "    [ -firstRun {FRUN} ] \\\n" \
  "    -runsPerBlock {NRB} \\\n" \
  "    -framesPerRun {NFR_PRE} {NFR_POS} \\\n" \
  "    { -stimulusChannel {ST_CHNAME[r]} {ST_TMAX[r]} {ST_OUTVAL[r]} {ST_RUNTYPE[r]} }.. \\\n" \
  "    -fixationChannel {FX_CHNAME} {FX_TMAX} {FX_OUTVAL} \\\n" \
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
  "  The program reads from standard input a 128-electrode EEG data file" \
  " from the (Bio/NonBio)x(Imb/Quiet) experiment by Ghislain et al [Gh2013b], and" \
  " extracts from it the segments that correspond to individual experimental" \
  " runs, as separate files.  It also removes unused trigger channels" \
  " and condenses the used ones to two new `fixation' and `stimulus' marker channels.\n" \
  "\n" \
  "INPUT FORMAT\n" \
  "  The input file should be a full-session EEG recording" \
  " of the 128-electrode equipment, converted to \".txt\" format.\n" \
  "\n" \
  "  The program assumes that the input file contains a certain" \
  " number {NT} of data frames, each with the same number {NC} of channel" \
  " values, including {NE} electrode readings and {NC-NE} non-electrode channels.\n" \
  "\n" \
  "  Some of the non-electrode channels are expected to be trigger channels, containing" \
  " pulses that signal the start and type of the `stimulus' phase of each" \
  " run.  Specifically, the stimulus phase begins when one of these" \
  " channels changes from 0 to a positive value.  These" \
  " channels are specified with the \"-stimulusChannel\" options.  An" \
  " additional channel, specified with the \"-fixationChannel\" option, is" \
  " assumed to contain pulses that signal the start of" \
  " the fixation phase that precedes the stimulus in each run.\n" \
  "\n" \
  "  The program assumes that each run includes specified frame counts {NFR_PRE} and {NFR_POS} before" \
  " and after the up-edge of the `stimulus' trigger.  Note that these counts may include" \
  " some frames at either end that are neither in the fixation phase nor in the stimulus" \
  " phase.  There is no explicit indication of the end of the stimulus phase, but it" \
  " is assumed to end with the start of the fixation phase of the next run, if any.\n" \
  "\n" \
  "  The runs in each file are conceptually grouped into `blocks' of" \
  " {NRB} consecutive runs, with variable delay between successive blocks.\n" \
  "\n" \
  "  The input file must have header records that specify, among other things, the number of" \
  " channels {NC} and their names {CHNAME[0]} {CHNAME[1]} ... {CHNAME[NC-1]}," \
  " the number of data frames {NT}, the sampling frequency {FSMP}, and" \
  " the number of electrode channels {NE}.  The first {NE} channels" \
  " in each data frame are assumed to be electrode" \
  " potentials, and the next one to be a `voltage reference' channel" \
  " called \"CZ\" that is always zero.  The channels with start-of-fixation" \
  " and start-of-stimulus pulses should be among the remaining {NC-NE-1} channels.\n" \
  "\n" \
  "OUTPUT FILES\n" \
  "  The output files are written with" \
  " names \"{OUT_DIR}/s{SSS}_r{BBB}{RR}.txt\" where {SSS}" \
  " is the three-digit subject ID number (obtained from the input file's header), {BBB} is" \
  " a three-digit sequential block number within the file, starting with {FBLOCK}, and {RR} is" \
  " a two-digit sequential run number within the block, starting with {FRUN}.  (See" \
  " the options \"-firstBlock\" and \"-firstRun\" below.)\n" \
  "\n" \
  "  Each output file contains the specified data frames surrounding a `stimulus'" \
  " pulse of the corresponding run.\n" \
  "\n" \
  "  Each output frame contains {NE+2} channels, comprising the same {NE} electrode" \
  " channels as in the input (excluding the reference electrode \"CZ\") followed by" \
  " two additional marker channels,\"FX\" and \"ST\", in that order.  The \"FX\" channel is set" \
  " to the specified value {FX_OUTVAL} throughout the fixation" \
  " phase of the run.  The \"ST\" channel is set to the" \
  " proper value {ST_OUTVAL[r]} (depending on stimulus type) throughout the stimulus phase.\n" \
  "\n" \
  "  Each output file begins with header records that specify some relevant" \
  " parameters such as number of frames and channels, sampling frequency, etc..\n" \
  "\n" \
  "OPTIONS\n" \
  "  -firstBlock {FBLOCK}\n" \
  "    This optional argument specifies the ID number of the first block" \
  " of runs in the input file, among all runs of the subject.  If" \
  " omitted, the program assumes {FRUN=1}.\n" \
  "\n" \
  "  -firstRun {FRUN}\n" \
  "    This optional argument specifies the ID number of the first run\n" \
  " in each block of runs.  If omitted, assumes {FRUN=1}.\n" \
  "\n" \
  "  -runsPerBlock {NRB}\n" \
  "    This mandatory argument specifies the number of runs in each" \
  " block of runs.\n" \
  "\n" \
  "  -framesPerRun {NFR_PRE} {NFR_POS}\n" \
  "    This mandatory argument specifies the number of frames to extract" \
  " for each run; specifically, {NFR_PRE} frames before the up-edge of" \
  " the start of `stimulus' pulse, and {NFR_POS} frames after it.\n" \
  "\n" \
  "  -stimulusChannel {ST_CHNAME[r]} {ST_TMAX[r]} {ST_OUTVAL[r]} {ST_RUNTYPE[r]}\n" \
  "    Each occurrence of this option keyword specifies the name {ST_CHNAME[r]} of" \
  " a trigger channel that may contain a start-of-stimulus pulse.  The" \
  " string {ST_CHNAME[r]} must match one of the" \
  " strings {CHNAME[NE+1..NC-1]} (excluding the first {NE+1} which are assumed" \
  " to be electrode voltages).\n" \
  "\n" \
  "    The stimulus that is marked with that channel will be" \
  " assumed to be of" \
  " type {ST_RUNTYPE[r]}.  The string {ST_RUNTYPE[r]} must start with" \
  " an ASCII letter and may contain only ASCII letters, digits" \
  " 0-9, periods '.' and underscores '_'.  This option must be specified at least once.\n" \
  "\n" \
  "    In the input file, the run's stimulus phase is assumed to begin with the" \
  " rise of channel {ST_CHNAME[r]} from 0 to any positive value and to end right" \
  " before the next start-of-fixation pulse, or after {ST_TMAX[r]} seconds, whichever comes first." \
  " In the output file, there will be a header line \"type = {ST_RUNTYPE[r]}\", and the new" \
  " marker channel \"ST\" will be set to {ST_OUTVAL} for the whole" \
  " duration of the stimulus phase.  \n" \
  "\n" \
  "  -fixationChannel {FX_CHNAME} {FX_TMAX} {FX_OUTVAL}\n" \
  "    This mandatory argument specifies the name {FX_CHNAME} of the" \
  " channel that contains the pulses signalling the start of" \
  " the `fixation' phase of each run.\n" \
  "\n" \
  "    The {FX_CHNAME} string must match one of the" \
  " strings {CHNAME[NE+1..NC-1]} (thus excluding the first {NE+1} which are assumed" \
  " to be electrode voltages), and must be distinct from all the `stimulus' channel" \
  " names {ST_CHNAME[r]}.\n" \
  "\n" \
  "    In the input file, the run's fixation phase is assumed to begin with the" \
  " rise of channel {FX_CHNAME} from 0 to any positive value and to end right" \
  " before the next start-of-stimulus pulse, or after {FX_TMAX} seconds," \
  " whichever comes first.  In the output" \
  " file, the new \"FX\" channel is set to {FX_OUTVAL} for the whole" \
  " duration of the fixation phase.\n" \
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
  "  Created 2013-10-03 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " nmeeg_split_e128_gh2013_C_COPYRIGHT ".\n" \
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
#include <neuromat_eeg_frame_buffer.h>

typedef struct nes_trigger_t
  { char *chName;   /* Name of channel with start of phase trigger pulses. */
    double tMax;    /* Max duration of phase in seconds. */
    double outVal;  /* New value for the marker channel in that phase. */
    char *runType;  /* Run type (NULL for fixation phase). */
  } nes_trigger_t;
  /* Specification of start-pulse channels for one run type. */
  
vec_typedef(nes_trigger_vec_t,nes_trigger_vec,nes_trigger_t);

typedef struct nes_options_t
  { 
    int firstBlock;        /* External number of first run block found in input file. */
    int firstRun;          /* External number of first run in each block. */
    int runsPerBlock;      /* Number of runs in each block. */
    int framesPerRun_pre;  /* Desired number of data frames before start-stimulus pulse. */
    int framesPerRun_pos;  /* Desired number of data frames after start-stimulus pulse. */
    char *outDir;          /* Directory for all output files. */
    nes_trigger_vec_t stimulusChannel; /* Trigger channel specs for each type of stimulus phase. */
    nes_trigger_t fixationChannel;     /* Trigger channel specs for the fixation phase. */
  } nes_options_t;
  /* Arguments from command line. */
  
#define nes_MAX_SUBJECT 999
  /* Max subject ID number. */

#define nes_MAX_BLOCK 99
  /* Max experimental blocks within a file. */

#define nes_MAX_RUNS_PER_BLOCK 100
  /* Max runs per block that the user may specify in the command line. */

#define nes_MAX_RUN (nes_MAX_RUNS_PER_BLOCK - 1)
  /* Max run ID number within an experimental block. */

#define nes_MAX_BLOCK_AND_RUN (nes_MAX_BLOCK*nes_MAX_RUNS_PER_BLOCK + nes_MAX_RUN)
  /* Max combined block ID and run ID number. */

neuromat_eeg_header_t *nes_read_header(FILE *rd, int *nlP);
  /* Reads from {rd} a neuroMat EEG dataset header, assuming the ".txt" 
    format, such as created by {neuromat_eeg_header_write}.
    
    If {nlP} is not null, {*nlP} is assumed to be the number of
    file lines (including headers and comments) previously read from
    {rd}; usually 0. The count is updated if not null. */
        
char **nes_define_output_chnames
  ( int nc_in, 
    char *chname_in[],
    int ne, 
    int np, 
    int ichs_out[],
    char *phnames_out[]
  );
  /* Returns a new list of names {chname_out[0..ne+np-1]} for 
    the output channels, comprising the {ne} input electrode channels
    plus {np} new phase marker channels.
    
    Entries {chname_out[0..ne-1]} are set to new clones of
    {chname_in[0..ne-1]}. Then each entry {chname_out[ichs_out[r]]} is
    set to a clone of {phnames_out[r]} for {r} in {0..np-1}. The indices
    {ichs_out[0..np-1]} must be a permutation of the new marker channel
    indices {ne..ne+np-1}. */

void nes_get_next_run
  ( neuromat_eeg_frame_buffer_t *buf,
    neuromat_eeg_frame_buffer_read_proc_t *read_frame,
    int it_start, 
    int ns, 
    int ichs[], 
    bool_t verbose,
    char *chname[], 
    double fsmp,
    int *it_fx_iniP,
    int *it_fx_finP, 
    int *it_st_iniP, 
    int *it_st_finP,
    int *s_stP
  );
  /* Looks for the next run in the input file, defined by a
    start-of-fixation pulse followed by a start-of-stimulus pulse,
    starting at or after frame {it_start}. Assumes that some frames are
    in {buf}, reads more with {read_frame} if needed. The frame
    {it_start} must not have been flushed from {buf} yet.

    Specifically, the procedure assumes that the valid start-of-stimulus pulses are
    in channels {ichs[0..ns-2]}, and that the start-of-fixation
    pulses are in channel {ic_fx_in=ichs[ns-1]} and. It looks for the first frame,
    beginning with {it_start}, where channel {ic_fx_in} is positive, and
    considers that the first frame of the start-of-fixation pulse. Then
    looks for the first frame after that where some channel {ic_st_in=ichs[s_st]}
    becomes positive, for some {s_st} in {0..ns-2}; and considers
    that the first frame of the start-of-stimulus pulse.

    If it succeeds, returns in {*it_fx_iniP,*it_fx_finP} the first and
    last frames of the start-of-fixation pulse, and in
    {*it_st_iniP,*it_st_finP} the first and last frames of the
    start-of-stimulus pulse. Also returns in {*s_stP} the the stimulus
    type in {0..ns-2}. Reports the pulse to
    {stderr}, assuming {chname[o..buf->nc-1]} are the channel names and
    {fsmp} is the sampling rate (in Hz).

    If there are no more runs in the input, returns {-1} in {*it_fx_iniP}.
    In that case {*it_fx_finP,*it_st_iniP,*it_st_finP,*s_stP} are undefined.

    Fails noisily if the next run is malformed in any way, e. g.
    if the first pulse found is a start-of-stimulus pulse, if the 
    following pulse is missing or another start-of-fixation pulse, 
    or if the two pulses overlap. */
    
double **nes_alloc_run(int nt, int nc);
  /* Allocate an EEG dataset {vru} with {nt} frames and {nc} channels,
    specifically {vru[0..nt-1][0..nc-1]}. */

void nes_extract_run
  ( neuromat_eeg_frame_buffer_t *buf,
    neuromat_eeg_frame_buffer_read_proc_t *read_frame, 
    int it_ini,
    int it_fin, 
    int nt,
    int nc,
    double **vru
  );
  /* Extracts a single experiment run from the EEG dataset buffer {buf},
    comprising frames {it_ini..it_fin}, and stores them in
    {vru[0..nt-1][0..nc-1]}.
    
    The procedure assumes that the values for each frame with index {it}
    in that interval are in {buf.val[itb][0..h.nc-1]} where
    {itb=(it%buf.size)}. It also requires {nt==it_fin-it_ini+1} and
    {nc==buf.nc}. The array {vru} must have been allocated by the
    client.
    
    Frame {it_ini} must not have been flushed from buffer {buf} yet. Reads
    more frames from {rd} if needed with the client-given procedure
    {read_frame}, which should return 1 if it succeeds, 0 if there
    are no more frames.  Fails if the file ends before frame {it_fin} */

void nes_get_phases
  ( int nt,
    int nc,
    double **vru,
    int ns,
    int ichs[], 
    int ic_fx_in,
    int nt_fx_max,
    int ic_st_in,
    int nt_st_max,
    int np,
    int kt_ini[],
    int kt_fin[]
  );
  /* Does some basic validity checking on the frames
    {vru[0..nt-1][0..nc-1]}.
    
    The procedure assumes that channels with indices {ichs[0..ns-1]}
    are stimulus trigger channels; that {ic_st_in}, in particular, is the stimulus
    trigger channel of this run; and that channel {ic_fx_in} is the 
    fixation trigger channel.  
    
    Within the frames {vru[0..nt-1]}, the trigger channels {ic_st_in} and
    {ic_fx_in} must define a non-empty fixation phase, immediately followed
    by a non-empty stimulus phase. Those two channels may not be
    simultaneously `on' (non-zero) in those frames,
    and all other trigger channels must be `off' (zero) in that interval.
    
    Specifically, the fixation phase frames {kt_fx_ini..kt_fx_fin} and
    the stimulus phase frames {kt_st_ini..kt_st_fin} are determined as
    follows:
    
      o Channel {ic_st_in} must be zero in frame {0}, and must have a
      single on-transition (from zero to non-zero) in {0..nt-1}.
      The first `on' frame is the start {kt_st_ini} of the stimulus
      phase.
    
      o Channel {ic_fx_in} may contain at most one on-transition before frame {kt_st_ini}.
      If there is such a transition, then the start {kt_fx_ini} of the fixation phase
      is the first `on' frame, otherwise it is frame {0}.
      
      o Channel {ic_fx_in} may contain at most one on-transition after frame {kt_st_ini}.
      If there is such a transition, then the end {kt_st_fin} of the stimulus phase
      is the last `off' frame; otherwise it is frame {nt-1}.
      
    In any case, the start {kt_fx_ini} of the fixation phase is
    incremented if needed so that the fixation it at most {nt_fx_max}
    frames long. Likewise the end {kt_st_fin} of the stimulus phase is
    decremented if needed so that the stimulus it at most {nt_st_max}
    frames long.
      
    The parameter {np} must be the number of expected phases, namely {2}.
    The frame index ranges for the two phases are returned in
    {kt_ini[0],kt_fin[0],kt_ini[1],kt_fin[1]}, respectively.
    Note that they relative to the start of the extracted run and not to
    the whole file. */

void nes_remap_channels
  ( int nt,
    int nc_in,
    double **vru,
    int ne,
    int np,
    int kt_ini[],
    int kt_fin[],
    int ichs_out[],
    double val_out[]
  );
  /* Remaps the input channels {vru[0..nt-1][0..nc_in-1]} to 
    to the output channels {vru[0..nt-1][0..nc_out-1]}.
    
    The procedure assumes that channels with indices {0..ne-1}
    are electrode measurements, that are to be preserved, 
    while channels in the range {ne..nc_out-1} are new 
    marker channels to be redefined. 
    
    Specifically, the procedure first sets
    {vru[0..nt-1][ne..nc_out-1]} to zeros;
    then, for each {r} in {0..np-1}, the procedure sets
    {vru[kt][ichs_out[r]]} to {val_out[r]}
    when {kt} is in {kt_ini[r]..kt_fin[r]}.
    
    The indices {ichs_out[0..np-1]} must be all distinct
    and in the range {ne..nc_out-1},
    and the ranges {kt_ini[r]..kt_fin[r]},
    if not empty, must be contained in {0..nt-1}. */

void nes_adjust_electrode_averages_in_run
  ( int nt,
    int nc,
    double **vru,
    int ic_fx_out,
    int ic_st_out,
    int ne
  );
  /* Shifts all electrode channels {0..ne-1} in the frames
    {vru[0..nt-1]} of a single run so that the `mean' value of
    each electrode is zero.  Returns the subtracted mean value.
    
    The averaging considers only samples in the fixation and stimulus phases
    (indicated by nonzero values in channels {ic_fx_out} and {ic_st_out},
    respectively) and omits apparent outliers. */

double nes_adjust_electrode_average
  ( int nt,
    int nc,
    double **vru,
    int ie,
    double wt[], /* (IN) Frame weights. */
    double pr[]  /* (WORK) Validity prob. */
  );
  /* Shifts the samples {vru[0..nt-1][ie]} so that the `mean' value of
    each electrode is zero.  In the mean computation, each frame {vru[kt]} 
    gets weight {wt[kt]}, and is further un-wheighted if its sample
    appears to be an outlier.  The scratch array {pr} must have {nt} elements. */
          
double nes_bayes(double v, double pri_gud, double avg_gud, double dev_gud, double avg_bad, double dev_bad);
  /* Uses Bayes's formula to compute the probability of {v} being a good sample,
    assuming good and bad samples have normal distributions with averages
    {avg_gud,avg_bad} and deviations {dev_gud,dev_bad}, respectively, 
    and that, a priori, the sample is good with probability {pri_gud}. */

void nes_compute_stats
  ( int nt, 
    int nc, 
    double **vru, 
    int ie, 
    double wt[], 
    double pr[], 
    bool_t good, 
    double *avgP, 
    double *devP, 
    double *priP
  );
  /* Estimates the average {*avgP} and deviation {*devP} of samples {vru[0..nt-1][ie]}.
    If {good} is TRUE, estimates the parameters of the `good' samples (inliers),
    if {good} is FALSE estimates those of the `bad' ones (outliers).
    Also returns in {*priP} (if not null) the estimated `a priori' probability of a 
    sample being good or bad, respectively.

    Assumes that {wt[kt]} is the relevance weight of sample {vru[kt][ie]},
    and {pr[kt]} is the current probability of that sample being a good sample. 
    If {pr} is null, assumes {pr[kt]==0.5} for all {kt}. */

void nes_compute_frame_weights_in_run(int nt, int nc, double **vru, int ic_fx_out, int ic_st_out, double *wt);
  /* Computes sample averaging weights for each frame in the run.
    Assumes that the frames are {vru[0..nt-1][0..nc-1]}, and that the
    fixation and stimulus phases are indicated by nonzero values in
    channels {ic_fx_out} and {ic_st_out}, respectively. Currently, the weight is
    1 inside the fixation and stimulus phases, 0 elsewhere. */

void nes_output_run
  ( int nt,
    int nc,
    double **vru,
    char *outDir,
    neuromat_eeg_header_t *h
  );
  /* Writes to disk one run extracted from the EEG dataset, comprising
    frames frames {vru[0..nt-1][0..nc-1]}.
    
    The output file will be called "{outDir}/s{SUBJ}_r{RUN}.txt", where
    {SUBJ} is the three-digit value of {h.orig.subject}, and {RUN} is
    the four-digit value of {h.orig.run}.
    
    The procedure sets {h.nt=nt} and {h.nc=nc}, and assumes that the client has
    defined all other fields of {h} properly, including {h.ne}, {h.type},
    {h.chname[0..nc-1]}, {h.orig.it_ini}, {h.orig.it_fin},
    {h.orig.run}, and {h.orig.subject}.
    
    The procedure then uses {neuromat_eeg_header_write} to write the
    non-null information from {h}, then calls {neuromat_eeg_data_write}
    to write the data. */

void nes_print_run_info(FILE *wr, char *class, neuromat_eeg_header_t *h);
  /* Prints to {wr} a one-line summary for one experiment run. Assumes that the run has the
    given {h.type} and spans from frame {h.orig.it_ini} to frame {h.orig.it_fin} of the
    original data file. The string {class} is printed in the summary.
    The sampling frequency {h.fsmp} is used to convert frame indices into
    time coordinates. */

void nes_print_run_phases
  ( FILE *wr, 
    int it_ini, 
    int it_fin, 
    int np,
    char *name[],
    int kt_ini[], 
    int kt_fin[],
    double fsmp
  );
  /* Prints to {wr} the start and stop of the {np} significant phases (such as fixation
    and stimulus) in a run.  Assumes that the extracted run starts with frame number
    {it_ini} and ends with frame {it_fin} of the original data file.
    Assumes that phase {i} of the run has codename {name[i]} and spans frames {kt_ini[i]..kt_fin[i]}
    both relative to the extracted run.  The sampling frequency
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
    bool_t debug = FALSE;
    
    /* Print some options: */
    fprintf(stderr, "first experimental block in input file is number %03d\n", o->firstBlock);
    fprintf(stderr, "first run of each block is number %02d\n", o->firstRun);
    fprintf(stderr, "expecting %d runs per block\n", o->runsPerBlock);
    fprintf(stderr, "extracting %d data frames before stimulus onset\n", o->framesPerRun_pre);
    fprintf(stderr, "extracting %d data frames after stimulus onset\n", o->framesPerRun_pos);
    
    int nl = 0; /* Number of lines read so far from {stdin}. */
    int nf = 0; /* Number of data frames read so far from {stdin}. */

    /* Read the file header: */
    neuromat_eeg_header_t *hin = nes_read_header(stdin, &nl);

    fprintf(stderr, "output files are named \"%s/s%03d_r{BBB}{NN}.txt\"\n", o->outDir, hin->orig->subject);

    int nc_in = hin->nc; /* Total count of data channels (incl. vref, triggers). */
    int ne = hin->ne; /* Count of electrode channels (excl. vref). */
    int nt = hin->nt; /* Count of data frames expected in input file. */
    
    /* Number of channels in output files: */
    int nc_out = ne; /* To begin with. */

    /* Define marker channels in output files, increment {nc_out}: */
    int ic_fx_out = nc_out;  nc_out++;
    int ic_st_out = nc_out;  nc_out++;
    assert(nc_out <= nc_in);
    
    /* Some data for the benefit of scripts that parse the {stderr} output: */
    fprintf(stderr, "subject = %03d\n", hin->orig->subject);
    fprintf(stderr, "nc_in = %d\n", nc_in);
    fprintf(stderr, "nc_out = %d\n", nc_out);
    fprintf(stderr, "ne = %d\n", ne);
    fprintf(stderr, "fsmp = %f\n", hin->fsmp);

    /* Convert trigger channel names to channel indices: */
    int ns = o->stimulusChannel.ne + 1; /* Count of start-of-phase channels (stimulus or fixation). */
    nes_trigger_t mk[ns];  /* {mk[0..ns-1]} are the specs for the {ns} trigger channels. */
    int ichs_in[ns];   /* {ichs_in[s]} is the channel index of start-of-stimulus marks for type {mk[s].rts}. */
    int s;
    for (s = 0; s < ns; s++)
      { char *phtitle = NULL;  /* Name of phase. */
        int ic_out = -1;
        if (s < ns-1)
          { mk[s] = o->stimulusChannel.e[s];
            asprintf(&phtitle, "stimulus of type %s:", mk[s].runType); 
            ic_out = ic_st_out;
          }
        else
          { mk[s] = o->fixationChannel;
            phtitle = txtcat("fixation:",""); 
            ic_out = ic_fx_out;
          }
        ichs_in[s] = neuromat_eeg_find_channel_by_name(mk[s].chName, ne+1, nc_in-1, hin->chname, TRUE);
        fprintf(stderr, "%-20s", phtitle);
        fprintf(stderr, "  start pulse in channel %d = %-4s", ichs_in[s], hin->chname[ichs_in[s]]);
        fprintf(stderr, "  will be mapped to channel %d value %+9.5f\n", ic_out, mk[s].outVal);
        free(phtitle);
      }
    
    int ic_fx_in = ichs_in[ns-1]; /* Index of channel where start-of-fixation appears. */
    int nt_fx_max = (int)imax(1, (int)ceil(o->fixationChannel.tMax * hin->fsmp));
    double val_fx_out = o->fixationChannel.outVal;
    
    /* Phase information: */
    int np = 2;    /* Number of phases per run. */
    int ichs_out[2] = { ic_fx_out, ic_st_out }; /* Marker channels in output files. */
    char *phnames_out[2] = { "FX", "ST" }; /* Channel names for fixation and stimulus marker channels. */

    /* Allocate the extracted run frames: */
    int nt_run = o->framesPerRun_pre + o->framesPerRun_pos; /* Frames per run. */
    double **vru = nes_alloc_run(nt_run, nc_in);
    /* Allocate the EEG data frame buffer: */
    int nt_buf = 2*nt_run; /* Frames to keep in buffer; should be enough. */
    neuromat_eeg_frame_buffer_t *buf = neuromat_eeg_frame_buffer_new(nt_buf, nc_in);

    /* Header for output files: */
    neuromat_eeg_header_t *hot = neuromat_eeg_header_new();
    neuromat_eeg_header_merge(hot, hin);
    hot->chname = nes_define_output_chnames(nc_in, hin->chname, ne, np, ichs_out, phnames_out);

    /* Run scanning data: */
    int idBlock = o->firstBlock; /* ID number of next block. */
    int idRun   = o->firstRun;   /* ID number of next run within block. */
    int nb = 0;                  /* Number of complete blocks found. */
    int nr = 0;                  /* Number of runs found. */
    int nrb = 0;                 /* Number of runs in current block. */
    int it_next = 0; /* Frame where to start looking for the next run. */
    bool_t verbose = TRUE;

    auto int read_frame(int nc_fr, double val_fr[]);
      /* Reads a frame from {rd} (in the plain text format)
        into {val_fr[0..nc_fr-1]}.  Expects {nc_fr==nc_in}. 
        Increments {nl} and {nf} as appropriate. */

    /* Main loop: Read from stdin splicing off the runs. */
    while (TRUE)
      { 
        /* Read another run from {stdin}: */
        if (debug) { fprintf(stderr, "trying to read block %d run %d\n", idBlock, idRun); }
        
        /* Find the start-of-fixation and start-of-stimulus pulses: */
        int it_fx_ini, it_fx_fin, it_st_ini, it_st_fin;
        int s_st;
        nes_get_next_run
          ( buf, read_frame, 
            it_next, ns, ichs_in, verbose, hin->chname, hin->fsmp, 
            &it_fx_ini, &it_fx_fin, &it_st_ini, &it_st_fin, &s_st
          );
        if (it_fx_ini < 0) { /* No more runs: */ break; }
          
        /* We got a new run: */
        demand((idBlock >= 0) & (idBlock <= nes_MAX_BLOCK), "invalid run block number");
        demand((idRun >= 0) & (idRun <= nes_MAX_RUN), "invalid run number");
        int idBlockRun = nes_MAX_RUNS_PER_BLOCK*idBlock + idRun; /* Run ID including block ID. */
        
        /* Some data for the benefit of scripts that parse the {stderr} output: */
        fprintf(stderr, "block = %03d\n", idBlock); /* For scripts that parse the log. */
        fprintf(stderr, "run = %02d\n", idRun); /* For scripts that parse the log. */
        fprintf(stderr, "blockrun = %05d\n", idBlockRun); /* For scripts that parse the log. */

        /* Dtermine the runtype and max stimulus frames from the trigger channel: */
        assert(s_st < ns-1);
        int ic_st_in = ichs_in[s_st];
        int nt_st_max = (int)imax(1, (int)ceil(mk[s_st].tMax * hin->fsmp));
        char *run_type = mk[s_st].runType;
        
        /* Check that the fixation phase was not too long: */
        demand(it_st_ini - it_fx_ini <= nt_fx_max, "fixation phase was too long");

        /* Decide the frame range {it_cp_ini..it_cp_fin} of frames to copy: */
        int it_cp_ini = it_st_ini - o->framesPerRun_pre;
        int it_cp_fin = it_st_ini + o->framesPerRun_pos - 1;

        /* Extract the run frames to {vru[0..nt-1][o..nc_in-1]}: */
        nes_extract_run(buf, read_frame, it_cp_ini, it_cp_fin, nt_run, nc_in, vru);

        /* Check trigger pulses, set marker channels to specified values: */
        int kt_ph_ini[2]; /* Initial frame of each phase (rel to start of run). */
        int kt_ph_fin[2]; /* Final frame of each phase (rel to start of run). */
        nes_get_phases
          ( nt_run, nc_in, vru, 
            ns, ichs_in, 
            ic_fx_in, nt_fx_max, 
            ic_st_in, nt_st_max,
            np,
            kt_ph_ini, kt_ph_fin 
          );
        nes_print_run_phases(stderr, it_cp_ini, it_cp_fin, np, phnames_out, kt_ph_ini, kt_ph_fin, hin->fsmp);

        /* Redefine the marker channels: */
        double val_st_out = mk[s_st].outVal;
        double val_out[2] = { val_fx_out, val_st_out };
        nes_remap_channels(nt_run, nc_in, vru, ne, np, kt_ph_ini, kt_ph_fin, ichs_out, val_out);
        
        /* Shift electrode averages to zero: */
        nes_adjust_electrode_averages_in_run(nt_run, nc_out, vru, ic_fx_out, ic_st_out, ne);

        /* Write run to disk: */
        hot->orig->run = idBlockRun;
        hot->type = run_type;
        hot->nc = nc_out;
        hot->orig->it_ini = it_cp_ini;
        hot->orig->it_fin = it_cp_fin;
        nes_output_run(nt_run, nc_out, vru, o->outDir, hot);
        
        /* Count one more run, prepare for next one: */
        nr++;
        nrb++;
        if (nrb >= o->runsPerBlock)
          { /* Start a new block: */
            idBlock++;
            idRun = o->firstRun;
            nb++;
            nrb = 0;
          }
        else
          { idRun++; }

        /* Assume that pulses do not overlap: */
        it_next = it_st_fin + 1;
      }
        
    /* Final reports: */
    fprintf(stderr, "read %d data frames in %d file lines\n", nf, nl);
    fprintf(stderr, "extracted total %d runs in %d complete blocks\n", nr, nb);
    if (nrb != 0) { fprintf(stderr, "!! warning: last block had only %d runs\n", nrb); }
    demand(nf == nt, "frame count in file does not match count in header");
    demand(nr > 0, "no pulses detected in input file");
    
    /* Release some storage: */ 
    int kt;
    for (kt = 0; kt < nt_run; kt++) { free(vru[kt]); } 
    free(vru);
    /* neuromat_eeg_frame_buffer_free(buf); */  /* Some bug here... */
    
    return 0;
    
    /* Internal implementations: */
    int read_frame(int nc_fr, double val_fr[])
      { assert(nc_fr == buf->nc);
        if (debug) { fprintf(stderr, "  trying to read frame %d at line %d\n", nf+1, nl+1); }
        int nrd = neuromat_eeg_frame_read(stdin, buf->nc, val_fr, &nl, &nf);
        assert((nrd == 0) || (nrd == 1));
        return nrd;
      }
  }

char **nes_define_output_chnames
  ( int nc_in, 
    char *chname_in[],
    int ne, 
    int np, 
    int ichs_out[],
    char *phnames_out[]
  )
  {
    int nc_out = ne + np; 
    char **chname_out = notnull(malloc(nc_out*sizeof(char*)), "no mem");
    int ic;
    /* Copy electrode names: */
    for (ic = 0; ic < ne; ic++) { chname_out[ic] = txtcat(chname_in[ic], ""); }
    /* Clear additional names: */
    for (ic = ne; ic < nc_out; ic++) { chname_out[ic] = NULL; }
    /* Set new channel names: */
    int r;
    for (r = 0; r < np; r++)
      { ic = ichs_out[r];
        demand((ic >= 0) && (ic < nc_out), "invalid output channel index");
        demand(chname_out[ic] == NULL, "duplicate output channel index");
        chname_out[ic] = txtcat(phnames_out[r], ""); 
      }
    /* Check if all are set: */
    for (ic = ne; ic < nc_out; ic++)
      { demand(chname_out[ic] != NULL, "undefined output channel"); }
    return chname_out;
  }

void nes_get_next_run
  ( neuromat_eeg_frame_buffer_t *buf,
    neuromat_eeg_frame_buffer_read_proc_t *read_frame,
    int it_start, 
    int ns, 
    int ichs[], 
    bool_t verbose,
    char *chname[], 
    double fsmp,
    int *it_fx_iniP,
    int *it_fx_finP, 
    int *it_st_iniP, 
    int *it_st_finP,
    int *s_stP
  )
  {
    /* In case of failure: */
    int it_fx_ini = -1;
    int it_fx_fin = -1;
    int it_st_ini = -1;
    int it_st_fin = -1;
    int s_st = -1;

    int s_fx = ns-1;  /* Phase type of the fixation phase: */
    
    /* Data for a start-of-phase pulse: */
    int it_pu_ini, it_pu_fin; /* First and last frame of pulse. */
    int s_pu;  /* Phase type in {0..ns-1}. */ 
    int ic_pu; /* Index of channel where pulse actually occurred. */
    
    /* Next frame to read: */
    int it_next = it_start;

    while (TRUE)
      { 
        /* Look for the next start-of-phase pulse in the fixation or stimulus trigger channels: */
        neuromat_eeg_frame_buffer_find_next_pulse
          ( buf, read_frame, it_next, 
            ns, ichs, 
            &it_pu_ini, &it_pu_fin, &s_pu
          );
        if (it_pu_ini < 0) 
          { /* End of file occurred before next start-of-phase pulse: */
            break;
          }
        assert((it_pu_ini >= it_next) && (it_pu_ini <= it_pu_fin));
        assert((s_pu >= 0) && (s_pu < ns));
        ic_pu = ichs[s_pu];
        if (verbose)
          { neuromat_eeg_report_pulse(stderr, " [fx] ", ic_pu, chname[ic_pu], it_pu_ini, it_pu_fin, fsmp, "\n"); }
        
        if (s_pu != s_fx)
          { fprintf(stderr, "!! warning: unexpected start-of-stimulus pulse in channel %d (ignored)\n", ic_pu);
            continue;
          }
        else
          { /* We got a start-of-fixation pulse: */
            it_fx_ini = it_pu_ini;
            it_fx_fin = it_pu_fin;
            break;
          }
      }
      
    if (it_fx_ini >= 0)
      { 
        /* Look for the start-of-stimulus pulse: */
        while (TRUE)
          { 
            /* Allow overlapping pulses for now. */
            it_next = it_fx_ini + 1; 

            /* Look for the next start-of-phase pulse in the fixation or stimulus trigger channels: */
            neuromat_eeg_frame_buffer_find_next_pulse
              ( buf, read_frame, it_next, 
                ns, ichs, 
                &it_pu_ini, &it_pu_fin, &s_pu);
            if (it_pu_ini < 0)
              { fprintf(stderr, "!! warning: file ended at frame %d after a fixation pulse (ignored)\n", buf->it_fin);
                /* Discard fixation pulse: */
                it_fx_ini = -1;
                it_fx_fin = -1;
                break;
              }
            assert((it_pu_ini >= it_next) && (it_pu_ini <= it_pu_fin));
            assert((s_pu >= 0) && (s_pu < ns));
            ic_pu = ichs[s_pu]; /* Channel where pulse actually occurred. */
            if (verbose)
              { neuromat_eeg_report_pulse(stderr, " [st] ", ic_pu, chname[ic_pu], it_pu_ini, it_pu_fin, fsmp, "\n"); }

            if (s_pu == s_fx)
              { fprintf(stderr, "!! warning: spurius start-of-fixation pulse at frame %d (ignored)\n", it_fx_ini);
                it_fx_ini = it_pu_ini;
                it_fx_fin = it_pu_fin;
              }
            else
              { /* We found the expected start-of-stimulus pulse: */
                it_st_ini = it_pu_ini;
                it_st_fin = it_pu_fin;
                s_st = s_pu;
                break;
              }
          }
      }
      
    if (it_fx_ini >= 0) 
      { assert(it_fx_fin >= it_fx_ini);
        assert(it_st_ini > it_fx_fin);
        assert(it_st_fin >= it_st_ini);
        assert((s_st >= 0) && (s_st < ns));
      }
    
    /* Return results: */
    (*it_fx_iniP) = it_fx_ini;
    (*it_fx_finP) = it_fx_fin;
    (*it_st_iniP) = it_st_ini;
    (*it_st_finP) = it_st_fin;
    (*s_stP) = s_st;
  }

double **nes_alloc_run(int nt, int nc)
  {
    double **vru = notnull(malloc(nt*sizeof(double*)), "no mem");
    int it;
    for (it = 0; it < nt; it++)
      { vru[it] = notnull(malloc(nc*sizeof(double)), "no mem"); }
    return vru;
  }

void nes_extract_run
  ( neuromat_eeg_frame_buffer_t *buf,
    neuromat_eeg_frame_buffer_read_proc_t *read_frame,
    int it_ini,
    int it_fin,
    int nt,
    int nc,
    double **vru
  )   
  {
    assert((nc > 0) && (nc == buf->nc));
    assert((nt > 0) && (nt == it_fin - it_ini + 1));
    
    /* Make sure frames {it_ini..it_fin} are in buffer: */
    int it;
    for (it = it_ini; it <= it_fin; it++)
      { int kt = it - it_ini; /* Index of frame {it} in {vru}. */
        assert(kt < nt);
        int itb = neuromat_eeg_frame_buffer_get_frame(buf, it, read_frame);
        demand(itb >= 0, "file ended with incomplete run");
        assert(itb < buf->size);
        int ic;
        for (ic = 0; ic < nc; ic++) { vru[kt][ic] = buf->val[itb][ic]; }
      }
  }
      
void nes_get_phases
  ( int nt,
    int nc,
    double **vru,
    int ns,
    int ichs[],
    int ic_fx_in,
    int nt_fx_max,
    int ic_st_in,
    int nt_st_max,
    int np,
    int kt_ini[],
    int kt_fin[]
  )      
  {
    /* Paranoia checks: */
    assert((ic_fx_in > 0) && (ic_fx_in < nc));
    assert(nt_fx_max > 0);
    assert((ic_st_in > 0) && (ic_st_in < nc) && (ic_st_in != ic_fx_in));
    assert(nt_st_max > 0);
    assert((ns > 0) && (ns <= nc-2)); /* At least one electrode and the fixation trigger. */
    
    /* Phase frame index ranges in {vru}, relative to first frame of run: */
    int kt_fx_ini = -1, kt_fx_fin = -1; /* Fixation phase. */
    int kt_st_ini = -1, kt_st_fin = -1; /* Stimlus phase. */
    
    /* Scan the frames.  The {state} variable has the following meanings:
      0 = before first frame of run.
      1 = after first frame, before any pulses.
      2 = during the first fixation pulse, before stimulus.
      3 = after the first fixation pulse, before stimulus.
      4 = during the stimulus pulse.
      5 = after the stimulus pulse, before the second fixation pulse.
      6 = during the second fixation pulse.
      7 = after the second fixation pulse.
    */
    int state = 0; /* Initial state. */ 
    int kt;
    for (kt = 0; kt < nt; kt++)
      { double *vt = vru[kt];
        demand(vt[ic_fx_in] >= 0, "fixation trigger channel is negative");
        /* Check that channels {ic_st_in,ic_fx_in} are consistent with {state}, update {state}: */
        bool_t on_fx = (vt[ic_fx_in] > 0); /* Is fixation trigger channel 'on'? */
        bool_t on_st = (vt[ic_st_in] > 0); /* Is this run's stimulus trigger channel 'on'? */
        demand(!(on_fx & on_st), "fixation and stimulus pulses overlap");
        switch(state)
          {
            case 0: /* This is the first frame of run: */
              demand(!on_st, "stimulus trigger is nonzero at start of run");
              /* There may be a fixation pulse: */
              if (on_fx)
                { kt_fx_ini = kt; state = 2; }
              else
                { state = 1; }
              break;

            case 1: /* Not first frame, no pulses yet: */
              /* There may be a pulse on either channel: */
              if (on_fx)
                { kt_fx_ini = kt; state = 2; }
              else if (on_st)
                { kt_st_ini = kt; state = 4;}
              else
                { state = 1; }
              break;

            case 2: /* Previous frame was in the first fixation pulse: */
              /* If the fixation pulse ends, the stimulus pulse may begin right away: */
              if (on_fx)
                { state = 2; }
              else if (on_st)
                { kt_st_ini = kt; state = 4;}
              else
                { state = 3; }
              break;

            case 3: /* Previous frame was between first fixation pulse and stimulus pulse: */
              /* The stimulus pulse must begin before the next fixation: */
              demand(!on_fx, "fixation pulse follows fixation pulse");
              if (on_st)
                { kt_st_ini = kt; state = 4; }
              else
                { state = 3; }
              break;

            case 4: /* Previous frame was inside stimulus pulse: */
              /* If the stimulus pulse ends, the fixation pulse may begin right away: */
              if (on_st)
                { state = 4; }
              else if (on_fx)
                { kt_st_fin = kt-1; state = 6;}
              else
                { state = 5; }
              break;

            case 5: /* Previous frame was between stimulus pulse and second fixation pulse: */
              demand(!on_st, "stimulus pulse follows stimulus pulse");
              if (on_fx)
                { kt_st_fin = kt-1; state = 6; }
              else
                { state = 5; }
              break;

            case 6: /* Previous frame was inside second fixation pulse: */
              /* Fixation pulse may end but stimulus may not begin: */
              demand(!on_st, "two stimulus pulses in run");
              if (on_fx)
                { state = 6; }
              else
                { state = 7; }
              break;

            case 7: /* Previous frame was after the second fixation pulse. */
              /* No more pulses allowed: */
              demand(!on_fx, "two fixation pulses after stimulus");
              demand(!on_st, "two stimulus pulses in short succession");
              break;

            default:
              assert(FALSE);
          }
        /* All other channels must be zero: */
        int r;
        for (r = 0; r < ns; r++)
          { int kc = ichs[r];
            demand(vt[kc] >= 0, "trigger channel is negative");
            if ((kc != ic_fx_in) && (kc != ic_st_in))
              { /* Stimulus trigger is not of this run, msut be zero: */
                demand(vt[kc] == 0.0, "spurious pulse in another stimulus channel");
              }
          }
      }
      
    /* We must have determined the start of the stimulus: */
    assert(kt_st_ini >= 0);
    kt_fx_fin = kt_st_ini - 1; /* By definition. */
    
    /* If there was no fixation pulse before stimulus, assume it was before first extracted frame: */
    if (kt_fx_ini < 0) { kt_fx_ini = 0; }
    
    /* If there was no fixation pulse after stimulus, assume it extends to end of extracted frame: */
    if (kt_st_fin < 0) { kt_st_fin = nt - 1; }
    
    /* Trim the two phases to the specified max length: */
    if (kt_fx_ini < kt_st_ini - nt_fx_max) { kt_fx_ini = kt_st_ini - nt_fx_max; }
    if (kt_st_fin > kt_st_ini + nt_st_max - 1) { kt_st_fin = kt_st_ini + nt_st_max - 1; }
    
    demand(np == 2, "incorrect number of phases");
    kt_ini[0] = kt_fx_ini;
    kt_fin[0] = kt_fx_fin;
    kt_ini[1] = kt_st_ini;
    kt_fin[1] = kt_st_fin;
  }

void nes_remap_channels
  ( int nt,
    int nc_in,
    double **vru,
    int ne,
    int np,
    int kt_ini[],
    int kt_fin[],
    int ichs_out[],
    double val_out[]
  )      
  {
    int nc_out = ne + np;

    /* Paranoia checks: */
    int r, ic, kt;
    for (r = 0; r < np; r++)
      { ic = ichs_out[r];
        assert((ic >= ne) && (ic < nc_out));
        if (kt_ini[r] <= kt_fin[r])
          { assert((kt_ini[r] >= 0) && (kt_fin[r] < nt)); }
        assert(! isnan(val_out[r]));
      }
    
    /* Ser the phase marker channels: */
    for (kt = 0; kt < nt; kt++)
      { double *vt = vru[kt];
        /* Fill channels {ne..ne+np-1} with {NAN}: */
        for (ic = ne; ic < nc_out; ic++) { vt[ic] = NAN; }
        /* Set marker channels to 0 or {val_out[r]} if within the respective phases: */
        for (r = 0; r < np; r++)
          { ic = ichs_out[r];
            if ((kt >= kt_ini[r]) && (kt <= kt_fin[r])) 
              { vt[ic] = val_out[r]; }
            else
              { vt[ic] = 0.0; }
          }
        /* Check for coverage: */
        for (ic = ne; ic < nc_out; ic++) { assert(! isnan(vt[ic])); }
      }
  }

void nes_adjust_electrode_averages_in_run
  ( int nt,
    int nc,
    double **vru,
    int ic_fx_out,
    int ic_st_out,
    int ne
  )     
  {
    demand((ic_fx_out >= ne) && (ic_fx_out < nc), "invalid fixation marker index");
    demand((ic_st_out >= ne) && (ic_st_out < nc), "invalid stimulus marker index");
    double *wt = notnull(malloc(nt*sizeof(double)), "no mem");  /* Importance weights. */
    double *pr = notnull(malloc(nt*sizeof(double)), "no mem");  /* Validity probs. */
    nes_compute_frame_weights_in_run(nt, nc, vru, ic_fx_out, ic_st_out, wt);
    int ie;
    for (ie = 0; ie < ne; ie++)
      { double avg = nes_adjust_electrode_average(nt, nc, vru, ie, wt, pr);
        fprintf(stderr, "    electrode %3d avg %+12.4f\n", ie, avg);
      }
    free(wt);
    free(pr);
  }
    
void nes_compute_frame_weights_in_run(int nt, int nc, double **vru, int ic_fx_out, int ic_st_out, double *wt)
  { 
    demand((ic_fx_out >= 0) && (ic_fx_out < nc), "invalid fixation marker index");
    demand((ic_st_out >= 0) && (ic_st_out < nc), "invalid stimulus marker index");
    int kt;
    /* Assign weight 1 inside fixation and stimulus phases: */
    for (kt = 0; kt < nt; kt++)
      { double *vt = vru[kt];
        wt[kt] = (double)((vt[ic_fx_out] != 0.0) | (vt[ic_st_out] != 0.0));
      }
    /* !!! Perhaps de-emphasize near edges? */
  }
   
double nes_adjust_electrode_average
  ( int nt,
    int nc,
    double **vru,
    int ie,
    double wt[],
    double pr[]
  )     
  {
    demand((ie >= 0) && (ie < nc), "invalid channel index");
    
    /* Paramters for the overall distribution: */
    double avg_all = NAN;   /* Assumed average of true samples. */
    double dev_all = NAN;   /* Assumed deviation of true samples. */
    nes_compute_stats(nt, nc, vru, ie, wt, NULL, TRUE,  &avg_all, &dev_all, NULL);
    assert(! isnan(avg_all));
    assert(! isnan(dev_all));
    
    /* Parameters of good samples (updated at each iteration): */
    double avg_gud = avg_all;      /* Assumed average of good samples. */
    double dev_gud = 0.5*dev_all;  /* Assumed deviation of good samples. */
    double pri_gud = 0.5;          /* A priori probability of sample being good. */
    /* Parameters of outlier distribution (to be adjusted): */
    double avg_bad = avg_all;     /* Assumed average of true samples. */
    double dev_bad = 4*dev_all;   /* Assumed deviation of true samples. */

    int max_iters = 5; /* Max iterations. */
    int iter;
    for (iter = 0; iter < max_iters; iter++)
      { bool_t debug = FALSE;
        if (debug) { fprintf(stderr, "    iteration %d\n", iter); }
        /* Cook the deviations if the distributions are equal: */
        bool_t same_avg = (fabs(avg_gud - avg_bad) < 0.05*dev_all);
        bool_t same_dev = (fabs(dev_gud - dev_bad) < 0.05*dev_all);
        if (same_avg && same_dev) { dev_gud = dev_gud/1.5; dev_bad = dev_bad*1.5; }
        /* Decide the probability {pr[kt]} of each sample {vru[kt][ie]} being an inlier: */
        int kt;
        for (kt = 0; kt < nt; kt++) { pr[kt] = nes_bayes(vru[kt][ie], pri_gud, avg_gud, dev_gud, avg_bad, dev_bad); }
        double pri_bad = NAN;
        /* Re-estimate the 'good' (inlier) and 'bad' (outlier) population parameters: */
        nes_compute_stats(nt, nc, vru, ie, wt, pr, TRUE,  &avg_gud, &dev_gud, &pri_gud);
        nes_compute_stats(nt, nc, vru, ie, wt, pr, FALSE, &avg_bad, &dev_bad, &pri_bad);
        if (debug) { fprintf(stderr, "      gud:  avg = %+9.4f  dev = %9.4f  pri = %9.4f\n", avg_gud, dev_gud, pri_gud); }
        if (debug) { fprintf(stderr, "      bad:  avg = %+9.4f  dev = %9.4f  pri = %9.4f\n", avg_bad, dev_bad, pri_bad); }
        assert(fabs(pri_bad + pri_gud - 1.0) < 0.0001);
        /* Re-estimate the 'bad' (outlier) population parameters: */
      }
      
    /* Subtract average from signal: */
    int kt;
    for (kt = 0; kt < nt; kt++)
      { double *vt = vru[kt];
        vt[ie] -= avg_gud;
      }
    return avg_gud;
  }
        
void nes_compute_stats
  ( int nt, 
    int nc, 
    double **vru, 
    int ie, 
    double wt[], 
    double pr[], 
    bool_t good, 
    double *avgP, 
    double *devP, 
    double *priP
  )
  {
    /* Compute the average and overall probability: */
    double sum_wpv = 0;
    double sum_wp = 0;
    double sum_w = 0;
    int kt;
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
    for (kt = 0; kt < nt; kt++) 
      { double dk = vru[kt][ie] - avg;
        double wk = wt[kt];
        double pk = (pr == NULL ? 0.5 : (good ? pr[kt] : 1.0 - pr[kt]));
        sum_wpd2 += wk*pk*dk*dk;
      }
    double dev = sqrt(sum_wpd2/sum_wp); /* Should compensate for bias... */
    /* Return: */
    assert(! isnan(avg));
    assert(! isnan(dev));
    assert(! isnan(pri));
    (*avgP) = avg;
    (*devP) = dev;
    if (priP != NULL) { (*priP) = pri; }
  }

double nes_bayes(double v, double pri_gud, double avg_gud, double dev_gud, double avg_bad, double dev_bad)
  {
    double dk_gud = (v - avg_gud)/dev_gud;
    double dk_bad = (v - avg_bad)/dev_bad;
    if (dk_gud > 6.0)
      { return 0.0; }
    else if (dk_bad > 6.0)
      { return 1.0; }
    else
      { double pri_bad = 1 - pri_gud;
        double P_gud = exp(-dk_gud*dk_gud/2)/dev_gud/sqrt(2*M_PI); /* Prob of sample, assuming good. */
        double P_bad = exp(-dk_bad*dk_bad/2)/dev_bad/sqrt(2*M_PI); /* Prob of sample, assuming bad. */ 
        double PP_gud = pri_gud*P_gud;
        double PP_bad = pri_bad*P_bad;
        return PP_gud/(PP_gud + PP_bad);
      }
  }

neuromat_eeg_header_t *nes_read_header(FILE *rd, int *nlP)   
  {
    neuromat_eeg_header_t *h = neuromat_eeg_header_read(stdin, 128, 500, nlP);

    fprintf(stderr, "original file name = \"%s\"\n", h->orig->file);
    fprintf(stderr, "subject ID = %03d\n", h->orig->subject);
    fprintf(stderr, "sampling frequency = %.10g\n", h->fsmp);

    fprintf(stderr, "expects %d data frames\n", h->nt);
    fprintf(stderr, "each input frame has %d channels =", h->nc);
    int ic;
    for (ic = 0; ic < h->nc; ic++) 
      { if ((ic % 10) == 0) { fprintf(stderr, "\n "); }
        fprintf(stderr, " %s", h->chname[ic]);
      }
    fprintf(stderr, "\n");
    fprintf(stderr, "channels 0..%d (%d channels) are electrodes\n", h->ne-1, h->ne);
    return h;
  }

void nes_print_run_phases
  ( FILE *wr, 
    int it_ini, 
    int it_fin, 
    int np,
    char *name[],
    int kt_ini[], 
    int kt_fin[],
    double fsmp
  )
  {
    demand((it_ini >= 0) && (it_ini <= it_fin), "invalid run frame range");
    int nt = it_fin - it_ini + 1;
    int ip;
    for(ip = 0; ip < np; ip++) 
      { demand((kt_ini[ip] >= 0) && (kt_ini[ip] <= kt_fin[ip]) && (kt_fin[ip] < nt), "invalid phase frame range");
        int nt_phase = kt_fin[ip] - kt_ini[ip] + 1;
        fprintf(wr, "   phase %2d = %s:  %8d frames ( %10.4f s )", ip, name[ip], nt_phase, nt_phase/fsmp);
        fprintf(wr, "  %8d .. %8d", kt_ini[ip], kt_fin[ip]);
        fprintf(wr, " ( %10.4f - %10.4f s )\n", (kt_ini[ip] - 0.5)/fsmp, (kt_fin[ip] + 0.5)/fsmp);
      }
    fprintf(wr, "\n");
  }

void nes_print_run_info(FILE *wr, char *class, neuromat_eeg_header_t *h)
  { 

    int it_ini = h->orig->it_ini;
    int it_fin = h->orig->it_fin; 
    double fsmp = h->fsmp;

    int nt = it_fin - it_ini + 1; /* Number of data frames in run. */
    double ts_ini_g = (it_ini - 0.5)/fsmp; /* Start time of run (seconds) from start of orig file. */
    double ts_fin_g = (it_fin + 0.5)/fsmp; /* End time of run (seconds) from start of orig file. */

    fprintf(wr, "subject %3d run %3d type %-8s (%s)", h->orig->subject, h->orig->run, h->type, class);
    fprintf(wr, " - %8d samples ( %10.4f s )", nt, nt/fsmp);
    fprintf(wr, " %8d .. %8d", it_ini, it_fin);
    fprintf(wr, " ( %10.4f _ %10.4f s )\n", ts_ini_g, ts_fin_g);
  }

void nes_output_run
  ( int nt,
    int nc,
    double **vru,
    char *outDir,
    neuromat_eeg_header_t *h
  )
  {
    /* Set {nc,nt} in header: */
    assert(nt > 0);
    assert(nc > 0);
    h->nt = nt;
    h->nc = nc;

    /* Check header data for this run: */
    assert((!isnan(h->fsmp)) && (h->fsmp == h->orig->fsmp)); /* Data is assumed to be unfiltered. */
    assert((h->ne > 0) && (h->ne <= h->nc-2)); /* There are at least two marker channels. */
    assert(h->type != NULL);
    assert(h->chname != NULL);
    assert(h->component == NULL);
    assert(h->kfmax == INT_MIN); /* The {kfmax} is irrelevant for time-domain data.  */
    assert(isnan(h->flo0));
    assert(isnan(h->flo1));
    assert(isnan(h->fhi1));
    assert(isnan(h->fhi0));
    assert(h->finvert == INT_MIN);
    
    /* Check header data re original file: */
    assert((h->orig->file != NULL) && (strlen(h->orig->file) > 0));
    assert(h->orig->nt > 0);
    assert((0 <= h->orig->it_ini) && (h->orig->it_ini <= h->orig->it_fin) && (h->orig->it_fin < h->orig->nt));
    assert((!isnan(h->orig->fsmp)) && (h->orig->fsmp > 0));
    assert((h->orig->subject > 0) && (h->orig->subject <= nes_MAX_SUBJECT));
    assert((h->orig->run > 0) && (h->orig->run <= nes_MAX_BLOCK_AND_RUN));
    
    nes_print_run_info(stderr, "extracted", h);
    
    char *fname = NULL;
    asprintf(&fname, "%s/s%03d_r%05d.txt", outDir, h->orig->subject, h->orig->run);
    FILE *wr = open_write(fname, TRUE);
    
    neuromat_eeg_header_write(wr, h);
    
    neuromat_eeg_data_write(wr, nt, nc, vru, 0, nt - 1, 1);
    fclose(wr);
    free(fname);
  }

#define nes_MAX_FRAMES_PER_RUN 72000
  /* Max frames in each output run that the user may specify
    in the command line. */
  
#define nes_MAX_PHASE_LENGTH (10.0)
  /* Max duration in seconds of a fixation or stimulus phase. */
  
#define nes_MIN_MARKER_VAL 0.01
#define nes_MAX_MARKER_VAL 1000.0
  /* Min and max output marker values that user may specify 
    in command line. */

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
    
    if (argparser_keyword_present(pp, "-firstBlock"))
      { o->firstBlock = (int)argparser_get_next_int(pp, 1, nes_MAX_BLOCK); }
    else
      { o->firstBlock = 1; }
      
    if (argparser_keyword_present(pp, "-firstRun"))
      { o->firstRun = (int)argparser_get_next_int(pp, 1, nes_MAX_RUN); }
    else
      { o->firstRun = 1; }
      
    /* Parse trigger channel name option and look it up in channel name list: */
    argparser_get_keyword(pp, "-runsPerBlock");
    o->runsPerBlock = (int)argparser_get_next_int(pp, 1, nes_MAX_RUNS_PER_BLOCK);

    argparser_get_keyword(pp, "-framesPerRun");
    o->framesPerRun_pre = (int)argparser_get_next_int(pp, 1, nes_MAX_FRAMES_PER_RUN);
    o->framesPerRun_pos = (int)argparser_get_next_int(pp, 1, nes_MAX_FRAMES_PER_RUN);
    
    argparser_get_keyword(pp, "-outDir");
    o->outDir = argparser_get_next_non_keyword(pp);
    
    /* Parse trigger list: */
    o->stimulusChannel = nes_trigger_vec_new(10);
    int ns = 0;
    while(argparser_keyword_present(pp, "-stimulusChannel"))
      { char *chName = argparser_get_next_non_keyword(pp);
        double tMax = argparser_get_next_double(pp, 1.0e-100, nes_MAX_PHASE_LENGTH);
        double outVal = argparser_get_next_double(pp, nes_MIN_MARKER_VAL, nes_MAX_MARKER_VAL);
        char *runType = argparser_get_next_non_keyword(pp);
        nes_trigger_vec_expand(&(o->stimulusChannel),ns);
        o->stimulusChannel.e[ns] = 
          (nes_trigger_t) { .chName = chName, .tMax = tMax, .outVal = outVal, .runType = runType };
        ns++;
      }
    nes_trigger_vec_trim(&(o->stimulusChannel),ns);
    if (ns == 0)
      { argparser_error(pp, "must specify at least one \"-stimulusChannel\" option"); }

    argparser_get_keyword(pp, "-fixationChannel");
    o->fixationChannel.chName = argparser_get_next_non_keyword(pp);
    o->fixationChannel.tMax = argparser_get_next_double(pp, 1.0e-100, nes_MAX_PHASE_LENGTH);
    o->fixationChannel.outVal = argparser_get_next_double(pp, nes_MIN_MARKER_VAL, nes_MAX_MARKER_VAL);
    o->fixationChannel.runType = NULL;
  
    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }

vec_typeimpl(nes_trigger_vec_t,nes_trigger_vec,nes_trigger_t);
