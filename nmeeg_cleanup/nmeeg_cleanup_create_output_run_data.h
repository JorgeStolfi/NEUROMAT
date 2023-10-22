#ifndef nmeeg_cleanup_create_output_run_data_H
#define nmeeg_cleanup_create_output_run_data_H

/* Assempbe the output run data for {nmeeg_cleanup}. */
/* Last edited on 2017-10-05 22:34:04 by jstolfi */

#define _GNU_SOURCE

#include <bool.h>

#include <nmeeg_cleanup_run_data.h>

#include <neuromat_eeg_header.h>

nmeeg_cleanup_run_data_t *nmeeg_cleanup_create_output_run_data
  ( nmeeg_cleanup_run_data_t *run_in,
    int ne,
    bool_t garb[]
  );
  /* Modifies the input EEG data {run_in} by seting garbage channels to zero,
    adding the channels "CZS", "BL0", "BL1",..., and "BLMK", and subtracting the
    new reference potential from all channels.
    
    A channel {ie} in {0..ne-1} is assumed to be garbage iff {garb[ie]}
    is true. The other information needed for this operation is obtained
    from the fields{.blc,.blmk,.vavg} of {run}. The boolean flags of
    {run.garb} must be a subset of {garb}. */

#endif
