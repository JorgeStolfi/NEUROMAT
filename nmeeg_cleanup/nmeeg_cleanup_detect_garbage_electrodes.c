/* See {nmeeg_cleanup_detect_garbage_electrodes.h} */
/* Last edited on 2017-10-17 17:59:05 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>

#include <neuromat_eeg_header.h>
#include <neuromat_eeg_channel_stats.h>

#include <nmeeg_cleanup_run_data.h>

#include <nmeeg_cleanup_detect_garbage_electrodes.h>

bool_t nmeeg_cleanup_detect_garbage_electrodes
  ( nmeeg_cleanup_run_data_t *run,
    int nk,             /* Number of marker channels that define relevant frames. */
    int ic_mk[],        /* Indices of those channels. */
    bool_t skipBlinks,  /* If true, ignores frames that are flagged as blinking episodes. */
    int nm,             /* Min number of correlated electrodes. */
    double corrThresh   /* Correlation threshold. */
  )
  {
    int nt = run->h->nt; /* Number of frames in run. */
    int nc = run->h->nc; /* Total number of channels in {run.val}. */
    int ne = run->h->ne; /* Number of electrode potentials in {run.val}. */
    
    /* Consistency of {ic_mk[]}: */
    for (int ik = 0; ik < nk; ik++)
      { int ic = ic_mk[ik];
        demand((ic >= ne) && (ic < nc), "invalid marker channel");
      }
      
    /* Define frame weights for average and covariance: */
    double wt[nt]; /* Frame weights. */
    int nf_ok = 0; /* Frames with nonzero weight. */
    for (int it = 0; it < nt; it++)
      { if (nk == 0)
          { /* In principle, include the frame: */
            wt[it] = 1.0;
          }
        else
          { /* In principle, include the frame if the given marker channels are nonzero: */
            wt[it] = 0.0;
            for (int ik = 0; ik < nk; ik++)
              { int ic = ic_mk[ik];
                if (run->val[it][ic] > 0.0) { wt[it] = 1.0; }
              }
          }
        /* Optionally exclude blink episodes: */
        if (skipBlinks && (run->blmk[it] > 0.0)) { wt[it] = 0.0; }
        /* Count relevat frames: */
        if (wt[it] != 0.0) { nf_ok++; }
      }
    demand(nf_ok >= 2, "too few frames for correlation analysis");
    
    /* Compute channel averages: */
    double vshift[ne];
    for (int ie = 0; ie < ne; ie++)
      { vshift[ie] = neuromat_eeg_channel_stats_avg(nt, nc, run->val, wt, ie); }

    /* Compute the electrode covariance matrix: */
    double *Cv = neuromat_eeg_channel_stats_covariance_matrix(nt, ne, run->val, vshift, wt);

    int nchanged = 0; /* Number of channels that changed status. */
    for (int ie = 0; ie < ne; ie++)
      { /* Compute the max {nmax} correlations between electrode {ie} and other electrodes: */
        double corr_max[nm];
        for (int im = 0; im < nm; im++) { corr_max[im] = -1.0; }
        for (int je = 0; je < ne; je++)
          { if (je != ie)
              { /* Compute correlation coeff {corr} between {ie} and {je}: */
                double Cii = Cv[ie*ne + ie];
                double Cij = Cv[ie*ne + je];
                double Cjj = Cv[je*ne + je];
                assert((Cii >= 0.0) && (Cjj >= 0.0));
                double Q = sqrt(Cii*Cjj + 1.0e-40);
                assert(fabs(Cij) <= Q);
                double corr = Cij/Q;
                /* Insert in queue of best correlations: */
                int km = nm;
                while ((km > 0) && (corr > corr_max[km-1])) 
                  { if (km < nm) { corr_max[km] = corr_max[km-1]; }
                    km--;
                  }
                if (km < nm) { corr_max[km] = corr; }
              }
          }
        /* Electrode {ie} is garbage if it has less than {nm} well-correlated neighbors: */
        bool_t garb_new = (corr_max[nm-1] < corrThresh);
        /* Update status and count changes: */
        if (garb_new != run->garb[ie]) 
          { fprintf(stderr, "electrode %3d garbage status changed", ie);
            fprintf(stderr, " from %c to %c\n", "FT"[run->garb[ie]], "FT"[garb_new]);
            nchanged++; 
          }
        run->garb[ie] = garb_new;
      }
    return (nchanged == 0);
  }
