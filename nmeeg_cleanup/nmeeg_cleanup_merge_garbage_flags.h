#ifndef nmeeg_cleanup_merge_garbage_flags_H
#define nmeeg_cleanup_merge_garbage_flags_H

/* Combining bad-electrode indicator vectors for {nmeeg_cleanup}. */
/* Last edited on 2017-10-05 22:29:34 by jstolfi */

#define _GNU_SOURCE

#include <bool.h>

#include <nmeeg_cleanup_run_data.h>

bool_t nmeeg_cleanup_merge_garbage_flags
  ( int nr, 
    nmeeg_cleanup_run_data_t **run, 
    int ne, 
    bool_t garb[]
  );
  /* Sets {garb[0..ne-1]} to the union of the indicator vectors
    {run[ir].garb[0..ne-1]}, for {ir} in {0..nr-1}.
    Returns {TRUE} iff the vector {garb} did not change. */

#endif
