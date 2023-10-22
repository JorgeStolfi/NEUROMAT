#ifndef nmeeg_cleanup_run_data_H
#define nmeeg_cleanup_run_data_H

/* Representation of a dataset for {nmeeg_cleanup}. */
/* Last edited on 2017-10-17 22:07:24 by jstolfi */

#define _GNU_SOURCE

#include <bool.h>

#include <neuromat_eeg_header.h>

typedef struct nmeeg_cleanup_run_data_t
  { char *runid;                  /* ID of run within subject/session, including block ID. */
    neuromat_eeg_header_t *h;     /* Header of dataset. */
    double **val;                 /* Data samples per frame and channel. */
    /* The following fields are used for input runs, during the cleanup: */
    bool_t *garb;                 /* TRUE for electrodes that seem trashed. */
    int nb;                       /* Number of components in blink pattern. */
    double **blc;                 /* Coefficient channels of blink pattern components. */
    double **blv;                 /* Combined blink component patterns. */
    double *blmk;                 /* Blinking marker channel for each frame in run. */
    double *vavg;                 /* Electrode average signal for re-referencing. */
  } nmeeg_cleanup_run_data_t;
  /* EEG data for one run, dirty or cleaned. 
    
    The sample array has the form {val[0..h.nt-1][0..h.nc-1]} where {val[it][ic]}
    is the sample in channel {ic} in the frame with index {it}.
    
    The fields {.garb,.blc,.blv,.blmk,.vavg} are used during and after the
    cleanup procedure.
    
    The booleans {.garb[0..h.nc-1]} tell which channels seem to have
    garbage values, e. g. electrodes with bad contact with the scalp.
    
    The reals {.blc[0..h.nt-1][0..nb-1]} are the coefficients of the
    {.nb} components of the blinking potential patterns present in each
    frame.  The linear combination of the patterns with those
    coeeficients are stored in {.blv[0..h.nt-1][0..h.ne-1]}.
    The reals {.blmk[0..nt-1]} are a new marker channel that
    identifies time intervals containing blinks.
    
    The reals {.vavg[0..nt-1]} are a weighted average of electrodes
    that is meant to become the new referene potential for the other
    electrodes. */

nmeeg_cleanup_run_data_t *nmeeg_cleanup_run_data_new(char *runid, neuromat_eeg_header_t *h);
  /* Allocates a new {nmeeg_cleanup_run_data_t} record. Set the blockrun
    ID {.runid} and dataset header {.h} fields from the given arguments.
    The {.nb} fiels is set to zero. Does not allocate the sample array
    {.val} or the work tables ({.garb,.blc}, etc) leaving those fields
    {NULL}. */
    
void nmeeg_cleanup_run_data_free(nmeeg_cleanup_run_data_t *run);
  /* Reclaims all storage in {run}, including sample and work tables. 
    Does not reclaim {run->h} or {run->runid}. */

void nmeeg_cleanup_run_data_alloc_work_tables(nmeeg_cleanup_run_data_t *run, int nb);
  /* Allocates the work tables {run.garb,run.blc}, etc., according to {run.h}.
    Sets the number of blink components {run.nb} to {nb} and allocates
    the {.blc} array with {nb} columns. */

void nmeeg_cleanup_run_data_free_work_tables(nmeeg_cleanup_run_data_t *run);
  /* Reclaims the work tables of {run}, if they exists. */

nmeeg_cleanup_run_data_t *nmeeg_cleanup_run_data_read(char *inPrefix, char *runid);
  /* Reads an EEG dataset from file "{inPrefix}_r{runid}.txt".
    The {runid} string is copied to a freshly allocated area. 
    Allocates the {.val} table but leaves the rest of the 
    work tables as {NULL}. */

void nmeeg_cleanup_run_data_write(char *outPrefix, char *runid, nmeeg_cleanup_run_data_t *run);
  /* Writes the signals {val[0..nt-1][0..nc-1]}, comprising {nc}
    channels sampled at {nt} times, to the file "{outPrefix}_r{runid}.txt".
    Writes the header {h} and then calls {neuromat_eeg_data_write} to
    write the data. */

void nmeeg_cleanup_run_data_print_blinks(FILE *wr, nmeeg_cleanup_run_data_t *run, bool_t secs);
  /* Prints to {wr} the frames {it} where {run.blmk[it]}
    is nonzero, as a comma-separated list of intervals. If {secs} is
    true prints the frame times instead of frame indices. */

#endif
