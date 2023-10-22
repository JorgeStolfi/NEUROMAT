#ifndef nmeeg_cleanup_detect_garbage_electrodes_H
#define nmeeg_cleanup_detect_garbage_electrodes_H

/* Detecting blinks and garbage channels for {nmeeg_cleanup}. */
/* Last edited on 2017-10-17 17:56:50 by jstolfi */

#define _GNU_SOURCE

#include <bool.h>

#include <nmeeg_cleanup_run_data.h>

bool_t nmeeg_cleanup_detect_garbage_electrodes
  ( nmeeg_cleanup_run_data_t *run,
    int nk,             /* Number of marker channels that define relevant frames. */
    int ic_mk[],        /* Indices of those channels. */
    bool_t skipBlinks,  /* If true, ignores frames that are flagged as blinking episodes. */
    int nm,             /* Min number of correlated electrodes. */
    double corrThresh   /* Correlation threshold. */
  );
  /* Performs garbage channel detection for the given {run}.
    
    The procedure evaluates the covariances of electrode channels
    in the run, and assumes that any channel {ie} with low correlation to
    other nearby electrodes is garbage.  In that case, the channel is flagged
    by setting {run.garb[ie]} to true. 
    
    In principle, the procedure considers only frames where at least one of the channels {ic_mk[0..nk-1]}
    is non-zero.  If {nc} is zero, in principle it considers all frames.
    In either case, however, if {skipBlinks} is true, the procedure ignores
    any frame {it} where {run.blmk[it]} is non-zero.
    
    The procedure returns {TRUE} iff no electrodes had their 
    `garbage' status changed by the call. */

#endif 
