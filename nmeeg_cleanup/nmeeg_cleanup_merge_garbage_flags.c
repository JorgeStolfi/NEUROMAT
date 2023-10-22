/* See {nmeeg_cleanup_merge_garbage_flags.h} */
/* Last edited on 2017-10-17 16:44:13 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>

#include <neuromat_eeg_header.h>
#include <neuromat_eeg_io.h>

#include <nmeeg_cleanup_run_data.h>

#include <nmeeg_cleanup_merge_garbage_flags.h>

bool_t nmeeg_cleanup_merge_garbage_flags
  ( int nr, 
    nmeeg_cleanup_run_data_t **run, 
    int ne, 
    bool_t garb[]
  )
  { 
    bool_t garb_next[ne];
    for (int ie = 0; ie < ne; ie++) { garb_next[ie] = FALSE; }
    for (int ir = 0;  ir < nr; ir++)
      { nmeeg_cleanup_run_data_t *runi = run[ir];
        assert(runi->h->ne == ne);
        for (int ie = 0; ie < ne; ie++) 
          { if (runi->garb[ie]) { garb_next[ie] = TRUE; } }
      }
    bool_t converged = TRUE;
    for (int ie = 0; ie < ne; ie++) 
       { if (garb[ie] != garb_next[ie]) { converged = FALSE; }
         garb[ie] = garb_next[ie];
       }
    return converged;
  }
