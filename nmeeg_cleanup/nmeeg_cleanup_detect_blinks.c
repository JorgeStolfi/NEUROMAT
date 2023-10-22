/* See {nmeeg_cleanup_detect_blinks.h} */
/* Last edited on 2017-10-18 06:21:07 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>

#include <bool.h>
#include <affirm.h>
#include <rmxn.h>

#include <neuromat_eeg_header.h>
#include <neuromat_eeg_io.h>
#include <neuromat_eeg_pca.h>

#include <nmeeg_cleanup_run_data.h>

#include <nmeeg_cleanup_detect_blinks.h>

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
  )
  {
    bool_t debug = TRUE;
    
    int nt = run->h->nt; /* Number of frames in run. */
    int ne = run->h->ne; /* Number of electrode potentisls. */
    demand((run->blc != NULL) && (run->blv != NULL), "needs aux tables for run");

    /* Compute the coefficients of the patterns in all frames, stores in {work1}: */
    neuromat_eeg_pca_fit_patterns(nt, ne, run->val, nb, blP, blQ, run->blc, work1, NULL);

    /* Compute the change {dc[0..nt-1]} in coefficient norm and new blink markers: */
    double dc[nt]; /* Euclidean change in blink component per frame. */
    bool_t blmkRaw[nt]; /* Raw blink marker per frame. */
    for (int it = 0; it < nt; it++)
      { double *blc_old = run->blv[it]; /* Blink component of frame {it} from previous iteration. */
        double *blc_new = work1[it]; /* Newly computed blink component of frame {it}. */
        if (debug)
          { fprintf(stderr, "    %05d blc = ", it);
            for (int ib = 0; ib < nb; ib++) 
              { fprintf(stderr, " %10.6f", run->blc[it][ib]); }
            fprintf(stderr, "\n");
          }
        blmkRaw[it] = (rn_norm(ne, blc_new) > blThresh);
        dc[it] = rn_dist(ne, blc_old, blc_new);
        rn_copy(ne, blc_new, blc_old);
      }
      
    /* Spread out blink markers: */
    double blmk_new[nt]; /* New blink barker channel. */
    for (int it = 0; it < nt; it++) { blmk_new[it] = 0.0; }
    int itIni = INT_MIN, itFin = INT_MIN; /* Range of frames to be marked as blinking. */
    for (int it = 0; it < nt; it++)
      { /* Now, if {itIni >= 0} all frames between {itIni} and {itFin} inclusive are to be marked. */
        if ((itIni >= 0) && (itFin + 1 < it - blmkRad))
          { /* The range {itIni..itFin} will not be extended further: */
            for (int jt = itIni; jt <= itFin; jt++) { blmk_new[jt] = blmkVal; }
            itIni = INT_MIN; itFin = INT_MIN;
          }
        if (blmkRaw[it])
          { if (itIni < 0) 
              { itIni = it - blmkRad;
                if (itIni < 0) { itIni = 0; }
              }
            itFin = it + blmkRad;
            if (itFin >= nt) { itFin = nt - 1; }
          }
      }
    if (itIni >= 0)
      { for (int jt = itIni; jt <= itFin; jt++) { blmk_new[jt] = blmkVal; } }
    int nmk_changed = 0; /* Number of frames that changed blink status. */
    for (int it = 0; it < nt; it++)
      { if ((blmk_new[it] != 0.0) != (run->blmk[it] != 0.0)) { nmk_changed++; }
        run->blmk[it] = blmk_new[it];
      }

    /* Check whether blink coeffs converged in non-blinking frames: */
    int nmk_changed_max = blmkRad; /* Tolerated change in blink status. */
    bool_t converged = (nmk_changed <= nmk_changed_max);
    for (int it = 0; (it < nt) && converged; it++)
      { if ((run->blmk[it] == 0.0) && (dc[it] > chThresh))
          { converged = FALSE; }
      }
    return converged;
  }
