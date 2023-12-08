/* See {nmeeg_cleanup_run_data.h} */
/* Last edited on 2023-12-05 23:38:55 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>

#include <bool.h>
#include <affirm.h>
#include <jsfile.h>
#include <jsmath.h>

#include <neuromat_eeg_header.h>
#include <neuromat_eeg_io.h>

#include <nmeeg_cleanup_run_data.h>

nmeeg_cleanup_run_data_t *nmeeg_cleanup_run_data_new(char *runid, neuromat_eeg_header_t *h)
  { 
    nmeeg_cleanup_run_data_t *run = notnull(malloc(sizeof(nmeeg_cleanup_run_data_t)), "no mem");
    run->runid = runid;
    run->h = h;
    int nt = h->nt;
    int nc = h->nc;
    double **val = alloc_C_matrix(nt, nc);
    run->val = val;
    /* Allocate work tables later if needed: */
    run->garb = NULL;
    run->nb = 0;
    run->blc = NULL;
    run->blv = NULL;
    run->blmk = NULL;
    run->vavg = NULL;
    return run;
  }
  
void nmeeg_cleanup_run_data_free(nmeeg_cleanup_run_data_t *run)
  {
    nmeeg_cleanup_run_data_free_work_tables(run);
    if (run->val != NULL) { free_C_matrix(run->val,run->h->nt); }
    free(run);
  }

void nmeeg_cleanup_run_data_alloc_work_tables(nmeeg_cleanup_run_data_t *run, int nb)
  { int nt = run->h->nt;
    int nc = run->h->nc;
    int ne = run->h->ne;
    nmeeg_cleanup_run_data_free_work_tables(run); /* Just in case. */
    run->garb = notnull(malloc(nc*sizeof(bool_t)), "no mem");
    run->nb = nb;
    run->blc = alloc_C_matrix(nt, nb);
    run->blv = alloc_C_matrix(nt, ne);
    run->blmk = notnull(malloc(nt*sizeof(double)), "no mem");
    run->vavg = notnull(malloc(nt*sizeof(double)), "no mem");
  }
  
void nmeeg_cleanup_run_data_free_work_tables(nmeeg_cleanup_run_data_t *run)
  { if (run->garb != NULL) { free(run->garb); run->garb = NULL; }
    if (run->blc  != NULL) { free_C_matrix(run->blc, run->h->nt); run->blc = NULL; }
    if (run->blv  != NULL) { free_C_matrix(run->blv, run->h->nt); run->blv = NULL; }
    if (run->blmk != NULL) { free(run->blmk); run->blmk = NULL; }
    if (run->vavg != NULL) { free(run->vavg); run->vavg = NULL; }
  }
    
nmeeg_cleanup_run_data_t *nmeeg_cleanup_run_data_read(char *inPrefix, char *runid)
  {
    char *fname = NULL;
    asprintf(&fname, "%s_r%s.txt", inPrefix, runid);
    FILE *rd = open_read(fname, TRUE);

    
    /* Read the input file's header, define {nc,ne}: */
    int nl = 0;
    neuromat_eeg_header_t *h = neuromat_eeg_header_read(rd, 20, 600, &nl);
    int nt = h->nt; /* Number of frames. */ 
    demand(nt > 0, "empty run file or invalid {nt} in file");
    int nc = h->nc; /* Number of channels in input file.*/
    demand(nc > 0, "invalid {nc} in file");
    int ne = h->ne; /* Number of electrode signals in input file. */
    demand((0 < ne) && (ne <= nc), "invalid {ne} in file");
    demand (h->chname != NULL, "missing channel names in input file header");
    /* int ic; for (ic = 0; ic < nc; ic++) { fprintf(stderr, "channel %3d = %-4s\n", ic, h->chname[ic]); } */
        
    /* Read the EEG data, define {nt} (ignore header {nt} if any): */
    int nt_read = 0; /* Number of data frames in run. */
    double **val = neuromat_eeg_data_read(rd, 0, h->nt, nc, &nl, &nt_read);
    fprintf(stderr, "read %d data frames\n", nt_read);
    demand(nt_read == nt, "did not read all frames");
    fprintf(stderr, "run has %d channels including %d electrodes\n", nc, ne);
    demand(nt > 0, "no frames to process");
    
    nmeeg_cleanup_run_data_t *run = nmeeg_cleanup_run_data_new(runid, h);
    run->runid = runid;
    run->h = h;
    run->val = val;
    
    free(fname);
    return run;
  }
    
void nmeeg_cleanup_run_data_write(char *outPrefix, char *runid, nmeeg_cleanup_run_data_t *run)
  {
    char *fname = NULL;
    asprintf(&fname, "%s_r%s.txt", outPrefix, runid);
    FILE *wr = open_write(fname, TRUE);

    neuromat_eeg_header_write(wr, run->h);
    neuromat_eeg_data_write(wr, run->h->nt, run->h->nc, run->val, "%14.8e", 0, run->h->nt-1, 1);
    fflush(wr);
    free(fname); 
  }

void nmeeg_cleanup_run_data_print_blinks(FILE *wr, nmeeg_cleanup_run_data_t *run, bool_t secs)
  { 
    auto void print_interval(char *s, int i0, int i1);
      /* Prints the range {i0..i1} prefixed with {s}, as
        seconds or indices, accoring to {secs}. */
    
    int nt = run->h->nt;
    int it_ini = -1, it_fin = -1;  /* Current blink interval, or {-1}. */
    char *sep = ""; /* Separator to use before next interval ("" or ","). */
    for (int it = 0; it < nt; it++)
      { if (run->blmk[it] == 0.0)
          { /* Frame is marked `not blinking': */
            if (it_ini >= 0)
              { print_interval(sep, it_ini, it_fin); 
                it_ini = -1; it_fin = -1; sep = ",";
              }
          }
        else
          { /* Frame is marked `blinking': */
            if (it_ini < 0) { it_ini = it; }
            it_fin = it;
          }
      }
    if (it_ini >= 0) { print_interval(sep, it_ini, it_fin); }

    return;
              
    void print_interval(char *s, int i0, int i1)
      { fprintf(wr, "%s", sep);
        if (secs)
          { /* Print interval in seconds: */
            double t0 = (i0 + 0.5)/run->h->fsmp;
            double t1 = (i1 + 0.5)/run->h->fsmp;
            fprintf(wr, "%.4f-%.4f", t0, t1);
          }
        else
          { /* Print interval of frame indices: */
            fprintf(wr, "%d-%d", i0, i1);
          }
      }

  }
