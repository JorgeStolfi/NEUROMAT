/* See {nmeeg_cleanup_create_output_run_data.h} */
/* Last edited on 2017-10-17 23:01:55 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#include <bool.h>
#include <affirm.h>
#include <jsmath.h>

#include <neuromat_eeg_header.h>
#include <neuromat_eeg_io.h>

#include <nmeeg_cleanup_run_data.h>

#include <nmeeg_cleanup_create_output_run_data.h>

nmeeg_cleanup_run_data_t *nmeeg_cleanup_create_output_run_data
  ( nmeeg_cleanup_run_data_t *run_in,
    int ne,
    bool_t garb[]
  )
  {
    /* Grab parameters of input run: */
    int nt = run_in->h->nt;
    int nc_in = run_in->h->nc;
    int ne_in = run_in->h->ne;
    demand(ne == ne_in, "inconsistent electrode count");
    
    /* Number of blink coefficient channels: */
    int nb = run_in->nb;
    
    /* Initialize the new header: */
    neuromat_eeg_header_t *h_out = neuromat_eeg_header_new();
    neuromat_eeg_header_merge(h_out, run_in->h);
    
    /* Add new channels to header: */
    int nc_out = nc_in + 1 + nb + 1;  /* "CZS", "BL0", "BL1", ..., "BLMK". */
    int ne_out = ne_in + 1 + nb;
    int ie_CZS = neuromat_eeg_header_append_electrode_channel(h_out, "CZS");
    int ie_BL0 = -1;
    for (int ib = 0; ib < nb; ib++) 
      { char *ename = NULL;
        asprintf(&ename, "BL%d", ib);
        int ie_tmp = neuromat_eeg_header_append_electrode_channel(h_out, ename);
        if (ib == 0) { ie_BL0 = ie_tmp; }
      }
    int ic_BLMK = neuromat_eeg_header_append_marker_channel(h_out, "BLMK");
    assert(h_out->ne == ne_out);
    assert(h_out->nc == nc_out);
    assert(h_out->nt == nt);
    assert(ie_CZS != ie_BL0);
    assert(ie_BL0 + nb == ne_out);
    assert(ic_BLMK == nc_out-1);
    
    /* Crete and fill the new sample array: */
    double **val_out = alloc_C_matrix(nt, nc_out);
    for (int it = 0; it < nt; it++)
      { double *v_in = run_in->val[it];  /* Original samples. */
        double *bl_in = run_in->blv[it]; /* Fitted blink component. */
        double *v_out = val_out[it];
        for (int ic = 0; ic < nc_out; ic++)
          { /* Define {v_out[ic]} from appropriate source: */
            if (ic < ne_in)
              { /* Copy input electrodes, subtracting blink and clearing out the garbage ones: */
                int ie = ic;
                if (garb[ie])
                  { v_out[ic] = 0.0; }
                else
                  { v_out[ic] = v_in[ic] - bl_in[ic]; }
              }
            else if (ic == ie_CZS)
              { v_out[ic] = - run_in->vavg[ic]; }
            else if ((ic >= ie_BL0) && (ic < ie_BL0 + nb))
              { int ib = ic - ie_BL0;
                v_out[ic] = run_in->blc[it][ib];
              }
            else if ((ic >= ne_out) && (ic < ic_BLMK))
              { int ic_in = ic - (ne_out - ne_in);
                v_out[ic] = v_in[ic_in];
              }
            else if (ic == ic_BLMK)
              { v_out[ic] = run_in->blmk[it]; }
            else
              { assert(FALSE); }
          }
      }

    char *runid_out = txtcat(run_in->runid, "");
    nmeeg_cleanup_run_data_t *run_out = nmeeg_cleanup_run_data_new(runid_out, h_out);
    run_out->val = val_out;
    return run_out;
  }
