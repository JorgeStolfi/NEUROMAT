/* Last edited on 2013-05-30 21:36:56 by stolfilocal */
/* Splits a multiple-run experiment into individual runs. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <fget.h>
#include <affirm.h>
#include <jsfile.h>
#include <jsmath.h>

#include <neuromat.h>

void nse_write_eeg_run(char *origfile, char *prefix, int nr_skip, int ir, char *run_type, int nt_skip, int nt, int nc, char *ch_name[], int ne, double fsmp, double **val, int it_ini, int it_fin, int np, int itp1[], int itp0[]);
  /* Writes a segment of the eeg signal array {val[0..nt-1][0..nc-1]},
    comprising {nc} channels sampled at {nt} times, to the file
    "{prefix}-r{nr_skip+ir}.txt". The segment spans samples with indices
    {it_ini} to {it_fin}, inclusive. Assumes that the data frame
    {val[0]} actually has index {nt_skip} in the input file. Requires 
    {0<=it_ini<=it_fin< nt}.
    
    Writes header lines with the original file name {origfile}, the
    sampling frequency {fsmp}, the run number {ir}, the run type
    {run_type}, the counts of sample times {nt}, channels {nc}, and
    electrodes {ne} in the original file, the channel names
    {ch_name[0..nc-1]}, the segment's index range
    {it_ini+nt_skip..it_fin+nt_skip} in the input file and the corresponding
    times (in seconds) since the start of the original data file. Then
    calls {neuromat_write_eeg_signals} to write the data. */

int main(int argc, char **argv)
  {
    int minargs = 8;
    demand(argc >= minargs, "bad args");
    int nextarg = 1;
    char *origfile = argv[nextarg]; nextarg++;    /* Name of raw input file (for documentation). */
    char *prefix = argv[nextarg]; nextarg++;       /* Prefix for all output filenames. */
    int nt_skip = atoi(argv[nextarg]); nextarg++;  /* Number of data frames to skip in input file. */
    int nt = atoi(argv[nextarg]); nextarg++;       /* Number of data frames to read from input file. */
    int nt_run = atoi(argv[nextarg]); nextarg++;   /* Number of data frames per run, or 0 if free. */
    int nr_skip = atoi(argv[nextarg]); nextarg++;  /* External number of first run found in read data. */
    int nr = atoi(argv[nextarg]); nextarg++;       /* Number of exerimental runs in file. */
    int np = atoi(argv[nextarg]); nextarg++;       /* Number of trigger pulses per run. */
    
    /* Get run types: */
    char **run_type = notnull(malloc(nr*sizeof(char*)), "no mem");
    int ir;
    for (ir = 0; ir < nr; ir++) { run_type[ir] = argv[nextarg]; nextarg++; }
    
    demand(nextarg >= argc, "spurious arguments");
    
    fprintf(stderr, "original file name %s\n", origfile);
    fprintf(stderr, "output files are named \"%s_r{NNN}.txt\"\n", prefix);
    fprintf(stderr, "skipping %d input data frames\n", nt_skip);
    fprintf(stderr, "reading %d input data frames\n", nt);
    if (nt_run > 0) { fprintf(stderr, "writing %d data frames per run\n", nt_run); }
    fprintf(stderr, "first experimental run in segment is number %03d\n", nr_skip);
    fprintf(stderr, "expecting %d experimental runs in read data\n", nr);
    fprintf(stderr, "assuming %d trigger pulses per run\n", np);
    
    demand(nt_skip >= 0, "invalid {nt_skip}");
    demand(nt > 0, "invalid {nt}");
    demand(nt_run >= 0, "invalid {nt_run}");
    demand(nr_skip >= 0, "invalid {nr_skip}");
    demand(nr > 0, "invalid {nr}");
    demand(np > 0, "invalid {np}");
    
    /* Assumed sampling frequency (Hz): */
    double fsmp = 600.0;
        
    /* Get the channel names: */
    int ne = 20;
    int nc = ne + 1; /* Includes {ne} channels plus one trigger (?) channel. */
    char **ch_name = neuromat_e20_signal_names();

    /* Read the EEG data: */
    double **val = neuromat_read_eeg_signals(stdin, nt_skip, nt, nc);
    fget_skip_formatting_chars(stdin);
    if (fget_test_char(stdin, EOF))
      { fprintf(stderr, "reached the end of the input file\n"); }
    else
      { fprintf(stderr, "!! warning: end of input file was not reached\n"); }
    
    /* Locate the intial and final frames of each run (tight range): */
    int *itrun = notnull(malloc(nr*sizeof(int)), "no mem"); /* Index of first frame of run in {0..nt-1}. */
    int *jtrun = notnull(malloc(nr*sizeof(int)), "no mem"); /* Index of last frame of run in {0..nt-1}. */
    int *kp1 = notnull(malloc(nr*np*sizeof(int)), "no mem"); /* Indices of trigger-up frames in {0..nt-1}. */
    int *kp0 = notnull(malloc(nr*np*sizeof(int)), "no mem"); /* Indices of trigger-down frames in {0..nt-1}. */
    ir = 0; /* Index of next run. */
    int ic_trig = nc - 1; /* index of trigger channel. */
    int kp = 0; /* Number of trigger-up events seen since last run. */
    double vtr_prev = 0; /* Trigger value in previous frame. */
    int it;
    for (it = 0; it < nt; it++)
      { double vtr = val[it][ic_trig];
        if ((vtr_prev == 0) && (vtr > 0))
          { /* Trigger up-event: */
            assert(kp < np);
            kp1[ir*np + kp] = it;
            if (kp == 0)
              { /* Start of run: */
                demand(ir < nr, "too many trigger events");
                itrun[ir] = it;
              }
            kp++;
          }
        else if ((vtr_prev > 0) && (vtr == 0))
          { /* Trigger down-event: */
            demand(kp > 0, "incomplete trigger pulse at start of file");
            assert(kp <= np);
            assert(it > 0);
            kp0[ir*np + (kp-1)] = it-1;
            if (kp == np)
              { /* End of run: */
                assert(ir < nr);
                jtrun[ir] = it-1;
                kp = 0;
                ir++;
              }
          }
        vtr_prev = vtr;
      }
    demand(vtr_prev == 0, "incomplete trigger pulse at end of file");
    demand(kp == 0, "incomplete run at end of file");
    demand(ir == nr, "insufficient number of runs in file");
      
    /* Decide the max number of frames to take before and after any run: */
    double ideal_buffer_pre =  1.25; /* Ideal buffer stretch (seconds)  before each run. */
    double ideal_buffer_pos =  1.25; /* Ideal buffer stretch (seconds) after each run. */ 
    int nt_buffer_pre = floor(fsmp*ideal_buffer_pre + 0.5); /* Buffer frames before run (if {nt_run} not given). */
    int nt_buffer_pos = floor(fsmp*ideal_buffer_pos + 0.5); /* Buffer frames after run (in any case). */
    if (nt_run == 0) 
      { fprintf(stderr, "buffer frames taken before each run = %d (%.4fs)\n", nt_buffer_pre, nt_buffer_pre/fsmp);
        fprintf(stderr, "buffer frames taken after each run =  %d (%.4fs)\n", nt_buffer_pos, nt_buffer_pos/fsmp);
      }
    
    /* Write expanded runs: */
    for (ir = 0; ir < nr; ir++)
      { char *fname = NULL;
        asprintf(&fname, "%s_r%04d.txt", prefix, ir); 
        int *itp1 = &(kp1[ir*np]); /* Indices of trigger-up frames of this run (rel to start of {val}). */
        int *itp0 = &(kp0[ir*np]); /* Indices of trigger-down frames of this run (rel to start of {val}). */
        assert(itp1[0] == itrun[ir]);
        assert(itp0[np-1] == jtrun[ir]);
        /* Get the frame index range {it_ini..it_fin} to extract: */
        int it_ini, it_fin;
        if (nt_run == 0) 
          { /* Start of run is {nt_buffer_pre} before beginning of first trigger pulse: */
            it_ini = itrun[ir] - nt_buffer_pre;
            it_fin = jtrun[ir] + nt_buffer_pos;
          }
        else
          { /* Divide the extra frames alloted between the two ends: */
            int nt_buffer = nt_run - (jtrun[ir] - itrun[ir] + 1);
            demand(nt_buffer > 0, "cannot honor {nt_run}, actual run is too long");
            it_fin = jtrun[ir] + nt_buffer/2;
            it_ini = it_fin - nt_run + 1;
          }

        /* Check whether the buffer frames exist and are indeed idle: */
        demand(it_fin < (ir == nr-1 ? nt : itrun[ir+1]), "not enough idle frames after end of run");
        demand(it_ini >= (ir == 0 ? 0 : jtrun[ir-1]+1), "not enough idle frames before start of run");
        
        /* Write the selected frames: */        
        nse_write_eeg_run
          ( origfile, prefix, nr_skip, ir, run_type[ir], nt_skip, 
            nt, nc, ch_name, ne, fsmp, val, it_ini, it_fin, np, itp1, itp0
          );
        free(fname);
      }
      
    for (it = 0; it < nt; it++) { free(val[it]); } 
    free(val);
    free(itrun);
    free(jtrun);
    free(kp0);
    free(kp1);
    free(ch_name);
    return 0;
  }

void nse_write_eeg_run(char *origfile, char *prefix, int nr_skip, int ir, char *run_type, int nt_skip, int nt, int nc, char *ch_name[], int ne, double fsmp, double **val, int it_ini, int it_fin, int np, int itp1[], int itp0[])
  {
    int it_ini_g = nt_skip + it_ini; /* Index of first segment frame in original file. */
    int it_fin_g = nt_skip + it_fin; /* Index of last segment frame in original file. */

    double ts_ini_g = (it_ini_g - 0.5)/fsmp; /* Start time of run (seconds) from start of orig file. */
    double ts_fin_g = (it_fin_g + 0.5)/fsmp; /* End time of run (seconds) from start of orig file. */
    
    int ir_g = nr_skip + ir; /* External run number. */

    char *fname = NULL;
    asprintf(&fname, "%s_r%03d.txt", prefix, ir_g);
    FILE *wr = open_write(fname, TRUE);

    fprintf(wr, "origfile = %s\n", origfile);
    fprintf(wr, "run = %d\n", ir_g);
    fprintf(wr, "runtype = %s\n", run_type);
    fprintf(wr, "indexrange = %d %d\n", it_ini_g, it_fin_g);
    fprintf(wr, "timerange = %.4f %.4f\n", ts_ini_g, ts_fin_g);
    fprintf(wr, "nc = %d\n", nc);
    fprintf(wr, "channels =");
    int ic;
    for (ic = 0; ic < nc; ic++) { fprintf(wr, " %s", ch_name[ic]); }
    fprintf(wr, "\n");
    fprintf(wr, "ne = %d\n", ne);
    fprintf(wr, "fsmp = %14.8f\n", fsmp);
    
    /* Print run summary to {stderr}: */
    int nt_run = it_fin - it_ini + 1;
    fprintf(stderr, "run %3d %-8s", ir_g, run_type);
    fprintf(stderr, " - %8d samples ( %10.4f s )", nt_run, nt_run/fsmp);
    fprintf(stderr, " %8d .. %8d", it_ini, it_fin);
    fprintf(stderr, " ( %10.4f _ %10.4f s )\n", ts_ini_g, ts_fin_g);
    if (np > 0)
      { /* Print the trigger pulses in this run: */
        fprintf(stderr, "\n");
        int k;
        for (k = 0; k < np; k++) 
          { int nt_pulse = itp0[k] - itp1[k] + 1;
            fprintf(stderr, "   trigger pulse %2d:  %8d samples ( %10.4f s )", k, nt_pulse, nt_pulse/fsmp);
            fprintf(stderr, "  %8d .. %8d", itp1[k] - it_ini, itp0[k] - it_ini);
            fprintf(stderr, " ( %10.4f - %10.4f s )\n", (itp1[k] - 0.5 - it_ini)/fsmp, (itp0[k] + 0.5 - it_ini)/fsmp);
          }
        fprintf(stderr, "\n");
      }
    
    neuromat_write_eeg_signals(wr, nt, nc, val, it_ini, it_fin);
    fclose(wr);
    free(fname);
  }
