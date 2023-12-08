#define PROG_NAME "nmeeg_convert_raw"
#define PROG_DESC "Reads a raw multiple-run EEG data file, dumps it."
#define PROG_VERS "2013-09-29"

#define nmeeg_convert_raw_C_COPYRIGHT \
  "Copyright © 2013 by the State University of Campinas (UNICAMP)"

/* Last edited on 2023-12-05 23:37:35 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    [ -unit {UNIT} ] \\\n" \
  "    [ -skipFrames {SKIPFRAMES} ] \\\n" \
  "    [ -copyFrames {COPYFRAMES} ] \\\n" \
  "    -sourceFile {SOURCEFILE} \\\n" \
  "    -subject {SUBJECT} \\\n" \
  "    [ -run {RUN} ] \\\n" \
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
  "  The program reads from standard input a raw EEG data file" \
  " in \".raw\" format, and writes to standard output" \
  " an eqivalent file in the \".txt\" format used by other NeuroMat programs.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -unit {UNIT}\n" \
  "    This optional argument specifies the unit of measurement" \
  " to assume for integer input data samples.  It is irrelevant if the" \
  " samples in the file are floating-point values.  If omitted or" \
  " zero, the program assumes {UNIT=1}.\n" \
  "\n" \
  "  -skipFrames {SKIPFRAMES}\n" \
  "    This optional argument tells how many data frames to skip in the" \
  " input file.  If omitted or zero, no frames are skipped.\n" \
  "\n" \
  "  -copyFrames {COPYFRAMES}\n" \
  "    This optional argument specifies the number of data frames to read after" \
  " the skipped frames.  If omitted or negative, reading continues to the end of file.\n" \
  "\n" \
  "  -sourceFile {SOURCEFILE}\n" \
  "  -subject {SUBJECT}\n" \
  "  -run {RUN}\n" \
  "    These arguments specify the name of the raw input" \
  " file (mandatory), the subject ID number (mandatory, from 1) and the ID number of the" \
  " experimental run (optional, from 1).  If the output file contains multiple" \
  " runs, leave {RUN} unspecified or set it negative. These" \
  " parameters do not affect" \
  " the program, they are only stored in the output file header" \
  " for documentation.  In particular, the raw file is read from" \
  " standard input, not from the file {SOURCEFILE}.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  neuromat_eeg_plot_signals.sh(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2013-09-29 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " nmeeg_convert_raw_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <math.h>

#include <argparser.h>
#include <affirm.h>
#include <jsfile.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_io.h>
#include <neuromat_eeg_raw_io.h>
#include <neuromat_eeg_raw_header.h>

typedef struct nrr_options_t
  { 
    double unit;       /* Unit of measurement for integer samples. */
    int skipFrames;    /* Data frames to read. */
    int copyFrames;    /* Max data frames to read. */
    char *sourceFile;  /* name of source file, for documentation only. */
    int subject;       /* Index of subject. */
    int run;           /* Index of experimental run. */
  } nrr_options_t;
  /* Arguments from command line. */

nrr_options_t *nrr_parse_options(int argc, char **argv);
  /* Parses the command line arguments. */
  
int main(int argc, char **argv);
  /* Main prog. */

int main(int argc, char **argv)
  {
    nrr_options_t *o = nrr_parse_options(argc, argv);
    
    /* Read the raw file header {hr}: */
    neuromat_eeg_raw_header_t *hr = neuromat_eeg_raw_header_read(stdin);
    neuromat_eeg_raw_header_print(stderr, hr);

    /* Decide how many records to skip: */
    int skipGoal = o->skipFrames;
    demand(skipGoal <= hr->nt, "skipping past end of input file");
    assert(skipGoal >= 0);
    
    /* Decide how many records to copy: */
    int copyGoal = (o->copyFrames < 0 ? hr->nt - skipGoal : o->copyFrames);
    demand(skipGoal + copyGoal <= hr->nt, "not enough frames in input file");
    assert(copyGoal >= 0);
    
    /* Convert to a plain header and write it out: */
    char *sourceFile = o->sourceFile;
    int subject = o->subject;
    int run = (o->run <= 0 ? INT_MIN : o->run); /* {INT_MIN} means unknown/multiple run. */
    neuromat_eeg_header_t *ht = neuromat_eeg_raw_header_to_plain_header
      ( hr, sourceFile, skipGoal, copyGoal, subject, run );
    neuromat_eeg_header_write(stdout, ht);

    int nc = ht->nc; /* Number of channels (incl. ref/timing and event channels). */
    assert(nc == hr->nc + hr->nv);
    int nt = ht->nt; /* Number of frames expected in input file. */
    assert(nt == copyGoal);
    
    double val[nc];
    int it;
    int skipCount = 0; /* Read but skipped frames. */
    int copyCount = 0; /* Read and written frames. */
    for (it = 0; it < skipGoal + copyGoal; it++)
      { int nr = neuromat_eeg_raw_frame_read(stdin, hr->version, o->unit, nc, ht->chname, val); 
        demand(nr > 0, "** unexpected EOF");
        if (it >= skipGoal) 
          { neuromat_eeg_frame_write(stdout, nc, val, "%14.8e");
            copyCount++;
          }
        else
          { skipCount++; }
      }
    fflush(stdout);

    fprintf(stderr, "skipped %d frames, copied %d frames\n", skipCount, copyCount);
    
    return 0;
  }

#define nrr_MAX_FILE_FRAMES ((1 << 30) - 1)
#define nrr_MIN_UNIT (0.000001f)
#define nrr_MAX_UNIT (1000000.0f)
#define nrr_MAX_SUBJECT (999999)
#define nrr_MAX_RUN (999999)

nrr_options_t *nrr_parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    nrr_options_t *o = notnull(malloc(sizeof(nrr_options_t)), "no mem"); 
    
    /* Parse keyword parameters: */
    
    if (argparser_keyword_present(pp, "-unit"))
      { o->unit = argparser_get_next_double(pp, nrr_MIN_UNIT, nrr_MAX_UNIT); }
    else
      { o->unit = 1.0; }
   
    if (argparser_keyword_present(pp, "-skipFrames"))
      { o->skipFrames = (int)argparser_get_next_int(pp, 0, nrr_MAX_FILE_FRAMES); }
    else
      { o->skipFrames = 0; }
    
    if (argparser_keyword_present(pp, "-copyFrames"))
      { o->copyFrames = (int)argparser_get_next_int(pp, -1, nrr_MAX_FILE_FRAMES); }
    else
      { o->copyFrames = -1; }
       
    argparser_get_keyword(pp, "-sourceFile");
    o->sourceFile = argparser_get_next(pp);
    
    argparser_get_keyword(pp, "-subject");
    o->subject = (int)argparser_get_next_int(pp, 1, nrr_MAX_SUBJECT);
      
    if (argparser_keyword_present(pp, "-run"))
      { o->run = (int)argparser_get_next_int(pp, INT_MIN, nrr_MAX_RUN); }
    else
      { o->run = INT_MIN; }
      
    /* Parse positional parameters: */
    argparser_skip_parsed(pp);

    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }
