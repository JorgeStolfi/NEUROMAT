/* Last edited on 2023-10-21 21:56:02 by stolfi */
/* Computes the power spectra of electrode signals an EEG dataset. */

/* !!! Convert to argparser.h !!! */

/* !!! Use trend instead of apodizing !!! */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <fftw3.h>

#include <r2.h>
#include <affirm.h>
#include <jsfile.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_io.h>
#include <neuromat_eeg_header.h>
#include <neuromat_filter.h>
#include <neuromat_eeg_stats.h>

void nesp_write_eeg_spectra(char *prefix, int ne, int kfmax, double **pwr, int nt, double fsmp, neuromat_eeg_header_t *h);
  /* Writes the power spectra {pwr[0..ne-1][0..kfmax]} of {ne} electrode
    signals to the file "{prefix}_pwr.txt". 
    
    Writes the header {h} after fixing its fields {h->nt,h->nc} to {ne}
    and {h->nt,h->fsmp,h->kfmax} to {nt,fsmp,kfmax}. */

int main(int argc, char **argv)
  {
    int nargs = 3;
    demand(argc >= nargs, "bad args");
    char *prefix = argv[1];      /* Prefix for all output filenames. */
    int nskip = atoi(argv[2]);   /* Index of first data frame to use. */
    int nread = atoi(argv[3]);      /* Number of data frames to use. */
    
    fprintf(stderr, "output file is %s_pwr.txt\n", prefix);
    if (nskip > 0) { fprintf(stderr, "skipping %d frames\n", nskip); }
    if (nread > 0) { fprintf(stderr, "reading %d frames\n", nread); }

    /* Read the EEG header, provide defaults: */
    int nl = 0;
    neuromat_eeg_header_t *h = neuromat_eeg_header_read(stdin, 20, 600, &nl);
    
    /* Get important parameters: */
    int nc = h->nc; /* Number of channels per data frame (including triggers etc.). */
    int ne = h->ne; /* Assumes that the first {ne} channels are electrode potentials. */
    double fsmp = h->fsmp;
    char **chname = h->chname;
    fprintf(stderr, "input file has %d channels, including %d electrode potentials\n", nc, ne);
    fprintf(stderr, "input sampling frequency %.10g Hz\n", fsmp);

    /* Read the EEG data: */
    int nt = 0;
    double **val = neuromat_eeg_data_read(stdin, nskip, nread, nc, &nl, &nt);
    fprintf(stderr, "read %d data frames from %d lines\n", nt, nl);
    
    /* Adjust original dataframe range to the range actually read: */
    if (h->orig->it_ini < 0) { h->orig->it_ini = 0; } 
    h->orig->it_ini += nskip;
    h->orig->it_fin = h->orig->it_ini + nt - 1;

    /* Compute basic statistics of electrode signals: */
    double vavg[nc], vvar[nc], vmin[nc], vmax[nc];
    double vmin_min, vmax_max, vdev_max;  /* Maximum variance among electrode channels. */
    bool_t zeroMean = FALSE;
    neuromat_eeg_stats_per_channel(nt, nc, val, zeroMean, vavg, vvar, vmin, vmax, ne, &vmin_min, &vmax_max, &vdev_max);
    fprintf(stderr, "--- signal statistics ---\n");
    neuromat_eeg_stats_per_channel_print(stderr, nc, chname, vavg, vvar, vmin, vmax);
    fprintf(stderr, "range of electrode values = [ %8.5f _ %8.5f ]\n", vmin_min, vmax_max);
    fprintf(stderr, "maximum dev = %8.5f\n", vdev_max); 
    fprintf(stderr, "\n");
 
     /* Compute  the power spectra of individual electrodes in original file: */
    int kfmax = nt/2;
    bool_t verbose = FALSE;
    double **pwr = neuromat_filter_compute_spectra(nt, ne, val, kfmax, verbose);
    nesp_write_eeg_spectra(prefix, ne, kfmax, pwr, nt, fsmp, h);
    
    double rpmin = 1.0e-12; /* Remove spectrum entries less than this rel total power. */
    double kfkeep = 0;      /* ... but keep at least up to this freq index. */
    if (rpmin > 0)
      { /* Adjust {kfmax} to range actually present in spectrum: */
        int ie;
        while (kfmax > kfkeep)
          { double trp = 0; /* Tot relative power at freq {kfmax}. */
            for (ie = 0; ie < ne; ie++)
              { double pwki = pwr[ie][kfmax];
                assert(pwki >= 0);
                trp += pwki/(vvar[ie] + 1.0e-100);
              }
            if (trp > rpmin) { break; }
            kfmax--;
          }
        double fmax = kfmax*fsmp/nt;
        fprintf(stderr, "spectrum truncated at max frequency %d ( %.10f Hz)\n", kfmax, fmax);
      }
    nesp_write_eeg_spectra(prefix, ne, kfmax, pwr, nt, fsmp, h);

    int it, ie;
    for (it = 0; it < nt; it++) { free(val[it]); } 
    free(val);
    for (ie = 0; ie < ne; ie++) { free(pwr[ie]); } 
    free(pwr);

    return 0;
  }

void nesp_write_eeg_spectra(char *prefix, int ne, int kfmax, double **pwr, int nt, double fsmp, neuromat_eeg_header_t *h)
  {
    demand(2*kfmax <= nt, "{kfmax} too high");
    
    /* Adjust header for this file: */
    h->nc = ne; /* No trigger/marker channels in spectrum. */
    h->ne = ne;
    h->nt = nt;      /* Original frame count, for the sake of freq computation. */
    h->fsmp = fsmp;  /* Original sampling frequency for the sake of freq computation. */
    h->kfmax = kfmax;
    
    char *fname = NULL;
    asprintf(&fname, "%s_pwr.txt", prefix);
    FILE *wr = open_write(fname, TRUE);

    neuromat_eeg_header_write(wr, h);
    neuromat_filter_write_spectra(wr, nt, ne, kfmax, fsmp, pwr);
    fclose(wr);
    free(fname);
  }
