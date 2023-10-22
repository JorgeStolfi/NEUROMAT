#ifndef nmeeg_cleanup_detect_blinks_H
#define nmeeg_cleanup_detect_blinks_H

/* Detecting blinks and garbage channels for {nmeeg_cleanup}. */
/* Last edited on 2017-10-18 06:07:01 by jstolfi */

#define _GNU_SOURCE

#include <bool.h>

#include <nmeeg_cleanup_run_data.h>

bool_t nmeeg_cleanup_detect_blinks
  ( nmeeg_cleanup_run_data_t *run,
    int nb,             /* Number of components in the blink pattern. */
    double *blP,        /* Component matrix of the potential pattern of blink pulses. */
    double *blQ,        /* Fitting matrix (inverse of {blP*blP'}). */
    double blThresh,    /* Threshold for blinking maker signal. */
    double blmkVal,     /* "True" value to use in the {blmk} channel. */
    int blmkRad,        /* Expand the {blmk} merker by this many frames. */
    double chThresh,    /* Convergence threshold. */
    double **work1      /* (WORK) Work matrix, {nt} rows and {ne} columns. */
  );
  /* Detects blinks in the given {run}.  To be used as one step in the
    main loop of the iterative cleanup algorithm.
    
    Let {nt=run.h.nt} and {ne=run.h.ne}. 
    
    Assumes that the blinks are described by the {nb} rows of matrix {blP}.
    Specifically, component {ib} of the blink pattern has potential {blP[ib*ne + ie]}
    in the electrode channel {ie}, for {ib} in {0..nb-1} and {ie} in
    {0..ne-1}.
    
    Ideally the rows of {blP} are orthonormal, and zero on previously
    detected garbage channels (in particular on all channels {ie} such
    that {run.garb[ie]} is true; but the procedure does not depend on
    it.  The procedure requires the {nb} by {nb} matrix {blQ} which should be
    the inverse of {blP*blP'} where {'} denotes transposition.
    
    The procedure will decompose each frame vector
    {vi=run.val[it][0..ne-1]}, for {it} in {0..nt-1}, into a linear
    combination {wi} of the rows of {blP} and a vector {ui} that is
    orthogonal to all those rows.
    
    The procedure sets {run.blc[it][ib]} to the coeeficient of component
    {ib} in that linear combination {wi}.  Stores {wi} into {run.blv[it][0..ne-1]}. 
    
    The procedure then sets the synthetic marker channel {run.blmk[it]}
    to {blmkVal} iff the Euclidean norm of the {wi} component
    is greater than {blThresh}. It then spreads
    out the flagged frames by {blmkRad} frames in each direction.
    
    The procedure returns {TRUE} iff the blink component {wi} 
    changed by less than {chThresh} on all frames that 
    were finally marked as non-blinking.
    
    The array {work1} must have at least {run.h.nt} rows and {ne} columns.*/

#endif 
