/* See {nmeeg_cleanup_determine_blink_patterns.h} */
/* Last edited on 2017-10-18 05:52:11 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <rmxn.h>
#include <rn.h>

#include <neuromat_eeg_header.h>
#include <neuromat_eeg_channel_stats.h>
#include <neuromat_eeg_pca.h>
#include <nmeeg_cleanup_run_data.h>

#include <nmeeg_cleanup_determine_blink_patterns.h>

void nmeeg_cleanup_determine_blink_patterns
  ( int nr,             /* Number of runs. */
    int nc,             /* Number of channels in all runs. */
    int ne,             /* Number of electrode potentials in all runs. */
    nmeeg_cleanup_run_data_t *run[], /* The input data for each run. */
    bool_t onlyBlinks,  /* If true, considers only frames previously marked as blinks. */
    bool_t garb[],      /* Tells which potential channels are garbage. */
    int nb,             /* Number of components in the blink potential pattern. */
    double *blP         /* (OUT) Components of the blink electrode potential of pattern. */ 
  )
  {
    demand(nr >= 1, "needs at least one run");
    
    double *Cv = notnull(malloc(ne*ne*sizeof(double)), "no mem");
    for (int ije = 0; ije < ne*ne; ije++) { Cv[ije] = 0.0; }
    
    for (int ir = 0; ir < nr; ir++)
      { nmeeg_cleanup_run_data_t *runi = run[ir];
        int nt = runi->h->nt;
        assert(runi->h->nc == nc);
        assert(runi->h->ne == ne);

        /* Define frame weights for average and covariance: */
        double wt[nt]; /* Frame weights for this run. */
        for (int it = 0; it < nt; it++)
          { /* In principle, include the frame: */
            wt[it] = 1.0;
            /* Optionally exclude blink episodes: */
            if (onlyBlinks && (runi->blmk[it] == 0.0)) { wt[it] = 0.0; }
          }

        /* Compute channel averages: */
        double vshift[ne];
        for (int ie = 0; ie < ne; ie++)
          { vshift[ie] = neuromat_eeg_channel_stats_avg(nt, nc, runi->val, wt, ie); }

        /* Accumulate the electrode covariance matrix: */
        neuromat_eeg_channel_stats_accum_covariance_matrix(nt, ne, runi->val, vshift, wt, 1.0, 0.0, Cv);
        
      }
      
    /* Clear covariances of garbage channels: */
    int ne_ok = ne; /* Number of non-garbage channels: */
    for (int ie = 0; ie < ne; ie++)
      { if (garb[ie])
          { for (int je = 0; je < ne; je++)
              { Cv[ie*ne + je] = 0.0;
                Cv[je*ne + ie] = 0.0; 
              }
            ne_ok--;
          }
      }
    demand(ne_ok >= nb + 1, "too few non-garbage channels left");

    /* Find eigenvectors and their magnitudes: */
    fprintf(stderr, "computing eigenvectors...\n");
    int nv; /* Number of principal components (significant eigenvectors): */
    double *Ev = rmxn_alloc(ne,ne); /* The first {nv} rows are the selected principal components. */
    double emag[ne]; /* Square roots of eigenvalues. */
    double minMag = 5.0; /* Ignore components with Euclidean norm less than this (µV). */
    nv = neuromat_eeg_pca_eigen_decomp(ne, Cv, minMag, Ev, emag); /* {nv} is num of eigens actually computed. */
    assert(nv <= ne_ok);
    if (nv < ne_ok) 
      { fprintf(stderr, "found only %d out of %d significant eigencomponents\n", nv, ne_ok); }
    
    /* Take the {nb} largest components as the blink pattern: */
    demand(nb <= nv, "not enough eigenvectors found");
    for (int ib = 0; ib < nb; ib++) 
      { double *Evi = &(Ev[ib*ne]);
        double *blPi = &(blP[ib*ne]);
        rn_copy(ne, Evi, blPi);
      }
      
    free(Cv);
    free(Ev);
  }
