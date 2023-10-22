#ifndef nmeeg_cleanup_determine_blink_patterns_H
#define nmeeg_cleanup_determine_blink_patterns_H

/* Determines the electrode potential pattern of blinks for {nmeeg_cleanup}. */
/* Last edited on 2017-10-17 18:00:36 by jstolfi */

#define _GNU_SOURCE

#include <bool.h>

#include <nmeeg_cleanup_run_data.h>

void nmeeg_cleanup_determine_blink_patterns
  ( int nr,             /* Number of runs. */
    int nc,             /* Number of channels in all runs. */
    int ne,             /* Number of electrode potentials in all runs. */
    nmeeg_cleanup_run_data_t *run[], /* The input data for each run. */
    bool_t onlyBlinks,  /* If true, considers only frames previously marked as blinks. */
    bool_t garb[],      /* Tells which potential channels are garbage. */
    int nb,             /* Number of components in the blink potential pattern. */
    double *blP         /* (OUT) Components of the blink electrode potential of pattern. */ 
  );
  /* Computes the electrode potential pattern charateristic of a blink.
    
    On input, {run[0..nr-1]} must be the runs of a single subject and session.
    
    If {onlyBlinks} is true, the {.blmk} field of each run must be
    non-null.  The procedure considers a frame {it} only if {.blmk[it]} nonzero.
    Otherwise all frames in all runs are considered.
    
    If {garb[ie]} is true for some {ie} in {0..ne-1}, the procedure assumes
    that channel {ie} is meaningless garbage and will ignore it.
    
    On output, {blP[0..nb*ne-1]} will be a matrix with {ne} columns
    whose {nb} rows are the main {nb} components of the electrode
    pattern of the blink pulses. Specifically, {blP[ib*ne + ie]} will be
    the potential of electrode {ie} in pattern {ib}, for {ib} in
    {0..nb-1} and {ue} in {0..ne-1}. In particular, that value will be
    zero if {garb[ie]} is true. The rows of {blP} will be orthogonal
    over the non-garbage electrode channels. */

#endif
