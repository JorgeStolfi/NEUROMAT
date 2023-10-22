/* See {neuromat_eeg.h}. */
/* Last edited on 2023-10-22 04:52:45 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <fget.h>
#include <affirm.h>
#include <jsstring.h>
#include <neuromat_eeg.h>

void neuromat_eeg_get_channel_names(int ne, int nv, char **evname, int *ncP, char ***chnamesP)
  { int nc = 0; /* Number of electrodes plust trigger/reference channel. */
    char **chname = NULL; /* Electrode names and trigger/reference channel name. */
    if (ne == 20)
      { neuromat_eeg_get_20_channel_names(&nc, &chname); }
    else if (ne == 128)
      { neuromat_eeg_get_128_channel_names(&nc, &chname); }
    else
      { demand(FALSE, "invalid electrode count"); }
    if (nv > 0)
      { /* Append the event channel names: */
        chname = notnull(realloc(chname, (nc+nv)*sizeof(char *)), "no mem");
        int i;
        for (i = 0; i < nv; i++) 
          { chname[nc+i] = txtcat(evname[i], ""); }
        nc = nc + nv;
      }
    /* Return results: */
    (*ncP) = nc;
    (*chnameP) = chname;
  }

void neuromat_eeg_get_20_channel_names(int *ncP, char ***chnameP)
  { 
    int nc = 21;
    char **chname = notnull(malloc(nc*sizeof(char*)), "no mem");
    chname[ 0] = "F7";
    chname[ 1] = "T3";
    chname[ 2] = "T5";
    chname[ 3] = "Fp1";
    chname[ 4] = "F3";
    chname[ 5] = "C3";
    chname[ 6] = "P3";
    chname[ 7] = "O1";
    chname[ 8] = "F8";
    chname[ 9] = "T4";
    chname[10] = "T6";
    chname[11] = "Fp2";
    chname[12] = "F4";
    chname[13] = "C4";
    chname[14] = "P4";
    chname[15] = "O2";
    chname[16] = "Fz";
    chname[17] = "Cz";
    chname[18] = "Pz";
    chname[19] = "Oz"; 
    chname[20] = "TR"; /* Trigger signal. */
    /* Make all strings be newly allocated: */
    int ic;
    for (ic = 0; ic < nc; ic++) { chname[ic] = txtcat(chname[ic], ""); }
    /* Return results: */
    (*ncP) = nc;
    (*chnameP) = chname;
  }
  
void neuromat_eeg_get_128_channel_names(int *ncP, char ***chnameP)  
  {
    int ne = 128;
    int nc = ne + 1;
    char **chname = notnull(malloc(nc*sizeof(char *)), "no mem");
    int i;
    /* Electrode channels: */
    for (i = 0; i < ne; i++) 
      { char *name = NULL;
        asprintf(&name, "C%03d", i+1);
        chname[i] = name;
      }
    /* Reference electrode: */
    chname[ne] = txtcat("CZ", ""); /* To get a newly allocated copy. */
    /* Return results: */
    (*ncP) = nc;
    (*chnameP) = chname;
  }

int neuromat_eeg_find_channel_by_name(char *name, int ic_start, int ic_end, char *chname[], bool_t die)
  {
    int ict = ic_start;
    while ((ict <= ic_end) && (strcmp(name, chname[ict]) != 0)) { ict++; }
    if (ict > ic_end)
      { if (die)
          { fprintf(stderr, "looking for channel \"%s\" in channels [%d..%d] = {", name, ic_start, ic_end);
            int jc;
            for (jc = ic_start; jc <= ic_end; jc++) { fprintf(stderr, " %s", chname[jc]); }
            fprintf(stderr, " }\n");
            demand(FALSE, "no such trigger channel");
          }
        else
          { ict = -1; }
      }
    return ict;
  }
    
int neuromat_eeg_locate_pulses
  ( int nt,              /* Number of data frames. */
    int nc,              /* Number of data channels. */
    double **val,        /* The EEG dataset ({nt} frames with {nc} channels each). */
    int ict,             /* Index of trigger channel. */
    double vmin,         /* "Low" trigger channel value. */
    double vmax,         /* "High" trigger channel value. */
    int np_max,          /* Max expected number of pulses in file. */
    int it_pulse_ini[],  /* (OUT) Index of first frame in each pulse. */
    int it_pulse_fin[]   /* (OUT) Index of last frame in each pulse. */
  )
  {
    demand((ict) && (ict < nc), "invalid trigger channel index");
    double vtr_prev = 0; /* Trigger channel value in previous frame. */
    int np = 0; /* Number of pulses found. */
    int it;
    for (it = 0; it <= nt; it++)
      { double vtr = (it < nt ? val[it][ict] : 0.0);
        demand((vtr == vmin) || (vtr == vmax), "invalid trigger value");
        if ((vtr_prev == vmin) && (vtr == vmax))
          { /* Trigger up-event: */
            demand(it > 0, "incomplete trigger pulse at start of file");
            demand(np < np_max, "too many trigger pulses");
            it_pulse_ini[np] = it;
            np++;
          }
        else if ((vtr_prev == vmax) && (vtr == vmin))
          { /* Trigger down-event: */
            demand(it < nt, "incomplete trigger pulse at end of file");
            assert(np > 0);
            assert(it > 0);
            it_pulse_fin[np - 1] = it - 1;
          }
        vtr_prev = vtr;
      }
    return np;
  }

void neuromat_eeg_report_pulse(FILE *wr, char *pre, int ic, char *name, int it_ini, int it_fin, double fsmp, char *suf)
  {
    if (pre != NULL) { fprintf(wr, "%s", pre); }
    fprintf(wr, "pulse in channel %d = \"%s\"", ic, name);
    fprintf(wr, " spanning frames %d..%d", it_ini, it_fin);
    fprintf(wr, " (%.4f _ %.4f sec)", ((double)it_ini)/fsmp, ((double)it_fin)/fsmp);
    if (suf != NULL) { fprintf(wr, "%s", suf); }
    fflush(wr);
  }
