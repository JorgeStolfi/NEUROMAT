#define PROG_NAME "nmeeg_comp_analysis"
#define PROG_DESC "Performs principal component analysis on an EEG dataset."
#define PROG_VERS "2013-06-06"

#define nmeeg_comp_analysis_C_COPYRIGHT \
  "Copyright © 2013 by the State University of Campinas (UNICAMP)"
/* Last edited on 2023-12-05 23:38:11 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    { -pattern {PAT_NAME} {PAT_FILE} }.. \\\n" \
  "    [ -normalize ] \\\n" \
  "    [ -writeCoeffs {COEFF_FILE} ] \\\n" \
  "    [ -writeComp {COMP_NAME} {COMP_FILE} = {PAT_NAME[0]} {PAT_NAME[1]} .. {PAT_NAME[M-1]} ].. \\\n" \
  "    [ -delete {CHNAME[0]} {CHNAME[1]} .. {CHNAME[M-1]} ] \\\n" \
  "    " argparser_help_info_HELP " \\\n" \
  "    < {INFILE} \\\n" \
  "    > {OUTFILE}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program reads from standard input an EEG dataset" \
  " and writes to standard output a version where new channels are defined by" \
  " certain linear combinations of original ones.  This progam can be used for example" \
  " to extract principal components identified by {nmeeg_correl}.\n" \
  "\n" \
  "  More specifically, the program assumes that every data line (/data frame/) in the input file contains" \
  " simultaneous samples from some fixed number {NC} of data channels, of which" \
  " the first {NE} are electrode potentials, while the rest are phase indicators (/markers/)" \
  " or other unspecified signals.\n" \
  "\n" \
  "  The program also reads a certain number {NP} of /pattern files/ {PAT[0..NP-1]} that define" \
  " the new channels.  Each pattern file should be an EEG dataset file with a single data" \
  " frame, with the same {NE} electrode channels as the input file.\n" \
  "\n" \
  "  The program views the first {NE} channel values in each input data" \
  " frame {t} as a vector {V[t][*] = V[t][NE-1]} in Cartesian space of dimension {NE}; and the" \
  " first {NE} values in each pattern file {PAT_FILE[k]} as a similar vector {P[k][*]}.  It" \
  " will decompose each input vector {V[t][*]} into two components: a linear combination {VP[t} of" \
  " the pattern vectors\n" \
  "\n" \
  "   {VP[t][*] = C[t][0]*P[0][*] + C[t][1]*P[1][*] + ··· + C[t][NP-1]*P[NP-1][*]}\n" \
  "\n" \
  " where {C[0..NP-1]} are {NP} real coefficients; and a residual" \
  " vector {VR[t][*]} with {NE} elements, that is orthogonal" \
  " to all vectors {P[0..NP-1]}.   The program will then concatenate the" \
  " vectors {VR[t][0..NE-1]} and {C[t][0..NP-1}" \
  " into an output vector {VO} with {NE+NP} channels.\n" \
  "\n" \
  "  To this vector the program will finally append the original {NC-NE} marker" \
  " channels from the input frame unchangde.  The resulting frame is then written to the output.\n" \
  "\n" \
  "  Optionally, some of the channels may be deleted before each frame is written out.\n" \
  "\n" \
  "COMPONENT FILES\n" \
  "  The program can also write separate data files that have the same {NE} electrodes as" \
  " the input, but whose values are combinations of some subset of the patterns found in the" \
  " original dataset, as fitted in each frame frame.  Namely, each frame {VS} of such an output file" \
  " will be a linear combination of pattern frames\n" \
  "\n" \
  "   {VS[t][*] = C[t][K[0]]*P[K[0]][*] + C[t][K[1]]*P[K[1]][*] + ··· + C[t][K[M-1]]*P[K[M-1]][*]}\n" \
  "\n" \
  " where {K[0..M-1]} is some subset of the pattern indices {0..NP-1}.\n" \
  "\n" \
  "FILE FORMAT\n" \
  "  The input file and the pattern files should have one or more header lines that" \
  " specify (among other parameters) the" \
  " number of channels {NC}, the number of electrodes {NE}, and the number of data" \
  " frames {NT}.  All datasets must have the same {NE} parameter.  In the input file, after" \
  " the header lines there must be {NT} data lines, each containing {NC} floating-point" \
  " numbers.  Each pattern file must contain a single data frame.\n" \
  "\n" \
  "  Output files have the same format, only with a different set of channels.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -pattern {PAT_NAME} {PAT_FILE}\n" \
  "    Each occurrence of this optional argument specifies the name {PAT_NAME} of a" \
  " pattern and the name of the file {PAT_FILE} that contains it.  The {PAT_NAME} must start with" \
  " an ASCII letter and may contain only ASCII letters, digits" \
  " 0-9, periods '.' and underscores '_'.  The {PAT_NAME}s must be distinct among themselves" \
  " and from other channel names.\n" \
  "\n" \
  "  -normalize\n" \
  "    This optional argument specifies that the electrode values of input patterns must be scaled after" \
  " reading so that the Euclidean norm of the electrode values is 1.  If this parameter is omitted, the" \
  " pattern frames are used as read.\n" \
  "\n" \
  "  -writeComp {COMP_NAME} {COMP_FILE} = {PAT_NAME[0]} {PAT_NAME[1]} .. {PAT_NAME[M-1]}\n" \
  "    Each occurrence of this optional argument requests the output to a" \
  " file called {COMP_FILE} of an EEG dataset" \
  " consisting of a linear combination of patterns {PAT_NAME[0..M-1]}, with the respective" \
  " coefficients as fitted in each frame.  Each {PAT_NAME[K]} must be the name of a" \
  " pattern defined in a preceding \"-pattern\" argument.   The {COMP_NAME} is a" \
  " string that will be written as the \"component\" field of the output" \
  " file's header; it must start with" \
  " an ASCII letter and may contain only ASCII letters, digits" \
  " 0-9, periods '.' and underscores '_'.\n" \
  "\n" \
  "  -writeCoeffs {COEFF_FILE}\n" \
  "    This optional argument requests the output to a" \
  " file called {COEFF_FILE} with in the  EEG dataset format" \
  " whose channels are the coefficients {C[0..NP-1}, plus the marker channels from the input.\n" \
  "\n" \
  "  -delete {CHNAME[0]} {CHNAME[1]} .. {CHNAME[M-1]}\n" \
  "    This optional argument requests that channels {CHNAME[0]} {CHNAME[1]} .. {CHNAME[M-1]}" \
  " be deleted before each frame is written out.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  neuromat_eeg_plot_signals.sh(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2013-12-01 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2013-12-01 J. Stolfi: split off from {nmeeg_correl}, substantial rewrite.\n" \
  "  2021-08-29 J. Stolfi: added \"-writeCoeffs\" option.\n" \
  "  2021-09-02 J. Stolfi: removed \"-norm\" option.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " nmeeg_comp_analysis_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define stringify(x) strngf(x)
#define strngf(x) #x

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <vec.h>
#include <rmxn.h>
#include <rn.h>
#include <argparser.h>
#include <affirm.h>
#include <jsstring.h>
#include <jsfile.h>
#include <jsmath.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_pca.h>
#include <neuromat_eeg_io.h>
#include <neuromat_eeg_header.h>
#include <neuromat_eeg_channel_stats.h>

typedef struct pat_spec_t
  { char *patname;  /* Name of new channel defined by the pattern. */
    char *fname;    /* Input file name. */
  } pat_spec_t;
  
vec_typedef(pat_spec_vec_t,pat_spec_vec,pat_spec_t);

typedef struct wcomp_spec_t
  { char *compname;         /* Name of channels to be combined. */
    string_vec_t patnames;  /* Name of channels to be combined. */
    char *fname;            /* Output file name. */
  } wcomp_spec_t;
  
vec_typedef(wcomp_spec_vec_t,wcomp_spec_vec,wcomp_spec_t);

typedef struct nca_options_t
  { pat_spec_vec_t pattern;     /* Channels defined by patterns. */
    wcomp_spec_vec_t writeComp; /* Components to write separately. */
    string_vec_t delete;        /* Channels to be deleted. */
    char *writeCoeffs;          /* File name for the coefficiesnts output file, or {NULL}. */
    bool_t normalize;           /* TRUE to normalize the patterns. */
  } nca_options_t;
  /* Arguments from command line. */
  
void nca_read_dataset(FILE *rd, neuromat_eeg_header_t **hP, double ***valP);
  /* Reads an EEG dataset from file {rd}, returns its sample array in {*valP} and the 
    header record {*hP}.  
    
    The sample array has the form {val[0..h.nt-1][0..h.nc-1]} where {val[t][i]}
    is the sample in channel {i} in the frame with index {t}. */
        
void nca_read_patterns
  ( int32_t ne,
    char *chname[],
    int32_t np, 
    pat_spec_t pattern[],
    bool_t normalize, 
    double **PP,
    double **QP
  );
  /* Reads {np} patterns from files specified in {pattern[0..np-1]}.
    Pattern {ip} is read from file {pattern[ip].fname}, which must have
    a single frame with {ne} electrodes, whose names must be
    {chname[0..ne-1]}.
    
    The procedure allocates and computes the pattern matrices {P} and
    {Q}, respectively of size {np} by {ne} and {np} by {np}, storing
    their addresses into {*PP} and {*QP}. The sample valies of those
    {ne} electrodes are stored in row {ip} of matrix {P}. */
  
void nca_read_pattern(char *fname, int32_t ne, char *chname[], bool_t normalize,  double pat[]);
  /* Reads a single-frame EEG dataset from file {fname}, returns the electrode
    values in {pat[0..ne-1]}.  The dataset must have at least {ne} channels
    whose names must match {chname[0..ne-1]}. Other channels are ignored. */
        
void nca_print_channel_stats(FILE *wr, int32_t nt, int32_t nc, int32_t ne, double **val, char *chname[]);
  /* Gathers and prints the per-channel statistics of the data set {val[0..nt-1][0..nc-1]}.
    Assumes that the first {ne} channels are electrode voltages (real or synthetic). */

void nca_delete_channels(int32_t nt, int32_t *ncP, int32_t *neP, double **val, char *chname[], bool_t keep[]);
  /* Assumes that {val} is an array with {nt} frames, each with {nc} channels,
    of thich the first {ne} are electrodes (real or synthetic);
    where {nc} is the input value of {*ncP}, and {ne} is the 
    input value of {*neP}.  Also assumes that {chname[0..nc-1]}
    are the (unshared) names of those channels, and that {keep} is an array of {nc}
    booleans.  Deletes from {chname} and from every frame of {val} 
    those channels with indices {ic} such taht {keep[ic]} is false.
    Compacts the remaining channels, keeping their order.  Then 
    updates {*ncP} and {*neP}with the number of surviving channels
    and surviving electrodes, respectively. */

void nca_write_eeg_dataset(FILE *wr, int32_t nt, int32_t nc, double **val, neuromat_eeg_header_t *h);
  /* Writes the signals {val[0..nt-1][0..nc-1]}, comprising {nc}
    channels sampled at {nt} times, to the file "{prefix}_{tag}.txt".
    Writes the header {h} and then calls {neuromat_eeg_data_write} to
    write the data. */

nca_options_t *nca_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments. */

string_vec_t nca_parse_chnames(argparser_t *pp);
  /* Parses a list of channel names from the command line, stopping at the next 
     keyword (any argument that starts with "-"). */


int32_t main(int32_t argc, char **argv);
  /* Main prog. */

int32_t main(int32_t argc, char **argv)
  {
    nca_options_t *o = nca_parse_options(argc, argv);
    
    /* Read the dataset: */
    neuromat_eeg_header_t *h = NULL; /* Parameters from the input file header. */
    double **vin = NULL;
    nca_read_dataset(stdin, &h, &vin);
    int32_t nt = h->nt;
    int32_t nc_in = h->nc;
    int32_t ne_in = h->ne;
    int32_t nm = nc_in - ne_in; /* Number of marker channels. */
    char **chname_in = h->chname;
    char **chname_in_mk = &(chname_in[ne_in]); /* Section of {chname_in} with marker names. */
    
    fprintf(stderr, "input channel statistics\n");
    nca_print_channel_stats(stderr, nt, nc_in, ne_in, vin, chname_in);
    
    int32_t np = o->pattern.ne; /* Number of pattern frames: */
    int32_t nd = o->delete.ne;  /* Number of channels to be deleted. */
    
    demand(np > 0, "no pattens given, aborted");
    
    /* Allocate the output dataset with all channels (before channel deletion): */
    int32_t nc_ot = nc_in + np; /* Number of temporary channels. */
    int32_t ne_ot = ne_in + np; /* Number of temporary electrodes. */
    double **vot = alloc_C_matrix(nt, nc_ot);
    
    /* Allocate the array of pattern coefficients, with extra space for marker channels: */
    double **coef = alloc_C_matrix(nt, np + nm);
        
    /* Read the patterns into {P[0..np*ne_in-1]}: */
    double *P = NULL; /* Pattern matrix. */
    double *Q = NULL; /* Matrix {P*P1}^{-1}}. */
    nca_read_patterns(ne_in, chname_in, np, o->pattern.e, o->normalize, &P, &Q);

    /* Creates the list of output channel names (including channels to be deleted): */
    char *chname_ot[nc_ot];
    char **chname_ot_pat = &(chname_ot[ne_in]); /* Section of {chname_ot} with the pattern names. */
    char **chname_ot_mk = &(chname_ot[ne_in + np]); /* Section of {chname_ot} with the merker names. */
    for (int32_t ie = 0; ie < ne_in; ie++) { chname_ot[ie] = txtcat(chname_in[ie], ""); }
    fprintf(stderr, "adding pattern channels:");
    for (int32_t ip = 0; ip < np; ip++) 
      { char *name = txtcat(o->pattern.e[ip].patname, "");
        fprintf(stderr, " %s", name);
        chname_ot_pat[ip] = name;
        int32_t kc = neuromat_eeg_find_channel_by_name(name, 0, nc_in-1, chname_in, FALSE);
        demand(kc < 0, "collision with input channel name");
        kc = neuromat_eeg_find_channel_by_name(name, 0, ip-1, chname_ot_pat, FALSE);
        demand(kc < 0, "collision with previous pattern name");
      }
    fprintf(stderr, "\n");
    for (int32_t im = 0; im < nm; im++) { chname_ot_mk[im] = txtcat(chname_in_mk[im], ""); }
    fprintf(stderr, "will have %d temporary channels:", nc_ot);
    for (int32_t ic = 0; ic < nc_ot; ic++) { fprintf(stderr, " %s", chname_ot[ic]); }
    fprintf(stderr, "\n");

    /* Fit the patterns to the data.
      Put the coeffs into {coef[0..nt][0..np-1]},
      and residue into in {vot[0..nt][0..ne_in-1]}: */
    fprintf(stderr, "fitting patterns to dataset\n");
    neuromat_eeg_pca_fit_patterns(nt, ne_in, vin, np, P, Q, coef, NULL, vot);

    /* Copy the coefficients into {vot} at channels {ne_in..ne_in+np-1}: */
    fprintf(stderr, "creating coefficients \"dataset\"\n");
    for (int32_t it = 0; it < nt; it++)
      { for (int32_t ip = 0; ip < np; ip++)
          { vot[it][ne_in + ip] = coef[it][ip]; }
      }

    /* Append the marker/trigger (non-electrode) channels from {vin} to {vot}: */
    for (int32_t ic = ne_in; ic < nc_in; ic++) 
      { for (int32_t it = 0; it < nt; it++) 
          { vot[it][ic+np] = vin[it][ic]; }
      }

    if (o->writeCoeffs != NULL)
      { /* Write the coefficients as a dataset: */
        fprintf(stderr, "writing  the coefficients \"dataset\"\n");
        char *coeffnames[np+nm];
        for (int32_t ip = 0; ip < np; ip++) { coeffnames[ip] = o->pattern.e[ip].patname; }
        for (int32_t im = 0; im < nm; im++) { coeffnames[np + im] =  chname_in[ne_in + im]; }

        /* Temporary header hacks: */
        h->ne = np;
        h->nc = np + nm;
        h->chname = coeffnames;
        h->component = "AMP";
        double *rebase_wt = h->rebase_wt;
        h->rebase_wt = NULL;

        FILE *wr = open_write(o->writeCoeffs, TRUE);
        nca_write_eeg_dataset(wr, nt, np + nm, coef, h);
        fclose(wr);

        /* Restore header: */
        h->ne = ne_in;
        h->nc = nc_in;
        h->chname = chname_in;
        h->component = NULL;
        h->rebase_wt = rebase_wt;
      }

    /* Write the requested fitted pattern combinations: */
    /* Warning: reuses {vin}. */
    int32_t nw = o->writeComp.ne;
    if (nw > 0)
      { for (int32_t iw = 0; iw < nw; iw++)
          { wcomp_spec_t *Siw = &(o->writeComp.e[iw]);  /* Spec of one output pattern combination file. */
            fprintf(stderr, "writing component %s to file %s\n", Siw->compname, Siw->fname);
            int32_t mp = Siw->patnames.ne; /* Number of patterns to combine into this file. */
            int32_t ip[mp]; /* Indices of patterns to combine. */
            for (int32_t k = 0; k < mp; k++) 
              { ip[k] = neuromat_eeg_find_channel_by_name(Siw->patnames.e[k], 0, np-1, chname_ot_pat, TRUE); }
            neuromat_eeg_pca_combine_patterns(nt, ne_in, np, P, coef, mp, ip, vin);
            h->component = Siw->compname;
            FILE *wr = open_write(Siw->fname, TRUE);
            nca_write_eeg_dataset(wr, nt, nc_in, vin, h);
            fclose(wr);
          }
      }

    if (nd > 0)
      { /* Delete channels: */
        bool_t keep[nc_ot]; /* TRUE for channels to keep. */
        for (int32_t ic = 0; ic < nc_in; ic++) { keep[ic] = TRUE; }
        for (int32_t id = 0; id < nd; id++)
          { char *delname = o->delete.e[id];
            int32_t ick = neuromat_eeg_find_channel_by_name(delname, 0, nc_ot-1, chname_ot, TRUE);
            demand(keep[ick], "duplicate channel in delete list");
            keep[ick] = FALSE;
          }
        nca_delete_channels(nt, &nc_ot, &ne_ot, vot, chname_ot, keep);
      }

    fprintf(stderr, "output channel statistics\n");
    nca_print_channel_stats(stderr, nt, nc_ot, ne_ot, vot, chname_ot);
    
    /* Write the output EEG dataset: */
    h->nc = nc_ot;
    h->chname = chname_ot;
    h->ne = ne_ot;
    h->component = NULL;
    nca_write_eeg_dataset(stdout, nt, nc_ot, vot, h);

    /* Dont's bother to free storage: */
    return 0;
  }

void nca_delete_channels(int32_t nt, int32_t *ncP, int32_t *neP, double **val, char *chname[], bool_t keep[])
  {
    int32_t nc_in = (*ncP);
    int32_t ne_in = (*neP);
    int32_t ic, nc_ot, ne_ot;
    
    /* Delete the channels in all frames: */
    int32_t it;
    for (it = 0; it < nt; it++)
      { double *v = val[it];
        for (ic = 0, nc_ot = 0; ic < nc_in; ic++)
          { if (keep[ic]) { v[nc_ot] = v[ic]; nc_ot++; } }
      }
    
    /* Delete the channel names: */
    for (ic = 0, nc_ot = 0, ne_ot = 0; ic < nc_in; ic++)
      { if (keep[ic]) 
          { chname[nc_ot] = chname[ic]; 
            if (ic < ne_in) { ne_ot++; }
            nc_ot++;
          } 
        else
          { free(chname[ic]); }
      }

    /* Update channel and electrode counts: */
    (*ncP) = nc_ot;
    (*neP) = ne_ot;
  }

void nca_read_dataset(FILE *rd, neuromat_eeg_header_t **hP, double ***valP)
  {
    int32_t nt = 0; /* Number of data frames read from all datasets. */
    
    /* Read the input file's header, define {nc,ne}: */
    int32_t nl = 0;
    neuromat_eeg_header_t *h = neuromat_eeg_header_read(rd, 20, 600, &nl);
    int32_t nc = h->nc; /* Number of channels in input file.*/
    demand(nc > 0, "invalid {nc} in file");
    int32_t ne = h->ne; /* Number of electrode signals in input file. */
    demand((0 < ne) && (ne <= nc), "invalid {ne} in file");
    demand (h->chname != NULL, "missing channel names in input file header");
    /* int32_t ic; for (ic = 0; ic < nc; ic++) { fprintf(stderr, "channel %3d = %-4s\n", ic, h->chname[ic]); } */
        
    /* Read the EEG data, define {nt} (ignore header {nt} if any): */
    double **val = neuromat_eeg_data_read(rd, 0, h->nt, nc, &nl, &nt);

    fprintf(stderr, "read %d data frames\n", nt);
    fprintf(stderr, "input file has %d channels including %d electrodes\n", nc, ne);
    demand(nt > 0, "no frames to process");
    
    (*hP) = h;
    (*valP) = val;
  }
    
  
void nca_read_patterns
  ( int32_t ne,
    char *chname[],
    int32_t np, 
    pat_spec_t pattern[],
    bool_t normalize, 
    double **PP,
    double **QP
  )
  { 
    demand(np <= ne, "too many patterns");

    /* Read the patterns and store into {P}: */
    double *P = rmxn_alloc(np,ne); /* Each row is a pattern frame. */
    for (int32_t ip = 0; ip < np; ip++)
      { pat_spec_t *Sip = &(pattern[ip]);
        double *Pip = &(P[ip*ne]);
        fprintf(stderr, "reading pattern %s from file %s\n", Sip->patname, Sip->fname);
        nca_read_pattern(Sip->fname, ne, chname, normalize, Pip);
      }

    /* Compute the array {Q[0..np*np-1]}, the inverse of {P*P'}: */
    fprintf(stderr, "computing the {Q} matrix\n");
    double *Q = rmxn_alloc(np,np); /* The array {(P*P')^-1}. */
    neuromat_eeg_pca_compute_fitting_matrix(np,ne,P,Q);
    fprintf(stderr, "done.\n");
    
    (*PP) = P;
    (*QP) = Q;
  }

void nca_read_pattern(char *fname, int32_t ne, char *chname[], bool_t normalize, double pat[])
  {
    FILE *rd = open_read(fname, TRUE); 
    neuromat_eeg_header_t *h = NULL;
    double **val = NULL;
    nca_read_dataset(rd, &h, &val);
    demand(h->ne == ne, "electrodes counts in pattern and file do not match");
    demand(h->nt == 1, "pattern should have a single data frame");
    demand (h->chname != NULL, "missing channel names in pattern header");
    double sum2 = 0;
    int32_t ie;
    for (ie = 0; ie < ne; ie++) 
      { demand(strcmp(chname[ie],h->chname[ie])== 0, "channel names in pattern and file do not match"); 
        double vi = val[0][ie];
        pat[ie] = vi;
        sum2 += vi*vi;
      }
    double norm = sqrt(sum2);
    fprintf(stderr, "pattern norm = %12.6f\n", norm);
    demand(norm > 1.0e-5, "pattern norm too small\n");
    if (normalize)
      { fprintf(stderr, "rescaling to unit rms norm\n");
        rn_scale(ne, 1/norm, pat, pat);
        fprintf(stderr, "done.\n");
      }
    free(val);
    free(h);
    fclose(rd);
  }
    
void nca_print_channel_stats(FILE *wr, int32_t nt, int32_t nc, int32_t ne, double **val, char *chname[])
  {
    /* Compute the channel means and variances over the selected frames: */
    neuromat_eeg_channel_stats_t *st = neuromat_eeg_channel_stats_new(nc);
    neuromat_eeg_channel_stats_t *stg = neuromat_eeg_channel_stats_new(1);
    double eps = 0.01; /* Assumed uncertainty of measurement (µV). */
    neuromat_eeg_channel_stats_gather_all(nt, nc, val, NULL, eps, st, ne, stg);
    fprintf(wr, "  --- channel statistics ---\n");
    neuromat_eeg_channel_stats_print_all(wr, 2, nc, chname, FALSE, st, ne, stg);
    fprintf(wr, "\n");
    free(st);
    free(stg);
  }

void nca_write_eeg_dataset(FILE *wr, int32_t nt, int32_t nc, double **val, neuromat_eeg_header_t *h)
  {
    assert(h->nt == nt);
    assert(h->nc == nc);
    fprintf(stderr, "writing final file...\n");
    for (int32_t ic = 0; ic < nc; ic++)
      { fprintf(stderr, "channel [%3d] = \"%s\"\n", ic, h->chname[ic]); }
    neuromat_eeg_header_write(wr, h);
    neuromat_eeg_data_write(wr, nt, nc, val, "%14.8e", 0, nt-1, 1);
    fflush(wr);
  }

nca_options_t *nca_parse_options(int32_t argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    /* Allocate the command line argument record: */
    nca_options_t *o = notnull(malloc(sizeof(nca_options_t)), "no mem");

    /* Parse keyword parameters: */

    { o->pattern = pat_spec_vec_new(10);
      int32_t np = 0;
      while (argparser_keyword_present(pp, "-pattern"))
        { char *patname = argparser_get_next_non_keyword(pp);
          char *fname = argparser_get_next_non_keyword(pp);
          pat_spec_vec_expand(&(o->pattern), np);
          o->pattern.e[np] = (pat_spec_t){ .patname = patname, .fname = fname };
          np++;
        }
      pat_spec_vec_trim(&(o->pattern), np);
    }
    
    { o->writeComp = wcomp_spec_vec_new(10);
      int32_t nw = 0;
      while (argparser_keyword_present(pp, "-writeComp"))
        { char *compname = argparser_get_next_non_keyword(pp);
          char *fname = argparser_get_next_non_keyword(pp);
          argparser_get_keyword_next(pp, "=");
          string_vec_t patnames = nca_parse_chnames(pp);
          wcomp_spec_vec_expand(&(o->writeComp), nw);
          o->writeComp.e[nw] = (wcomp_spec_t){ .compname = compname, .fname = fname, .patnames = patnames };
          nw++;
        }
      wcomp_spec_vec_trim(&(o->writeComp), nw);
    }

    if (argparser_keyword_present(pp, "-writeCoeffs"))
      { o->writeCoeffs = argparser_get_next_non_keyword(pp); }
    else
      { o->writeCoeffs = NULL; }
    
    { o->delete = string_vec_new(10);
      int32_t nd = 0;
      while (argparser_keyword_present(pp, "-delete"))
        { string_vec_t chname = nca_parse_chnames(pp);
          int32_t id;
          for (id = 0; id < chname.ne; id++)
            { string_vec_expand(&(o->delete), nd);
              o->delete.e[nd] = chname.e[id];
              nd++;
            }
        }
      string_vec_trim(&(o->delete), nd);
    }
    
    o->normalize = argparser_keyword_present(pp, "-normalize");

    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }


string_vec_t nca_parse_chnames(argparser_t *pp)
  { 
    string_vec_t chname = string_vec_new(10);
    int32_t nc = 0;
    while (argparser_next_is_non_keyword(pp))
      { char *name = argparser_get_next_non_keyword(pp);
        string_vec_expand(&chname, nc);
        chname.e[nc] = name;
        nc++;
      }
    string_vec_trim(&chname, nc);
    return chname;
  }

vec_typeimpl(pat_spec_vec_t,pat_spec_vec,pat_spec_t);
vec_typeimpl(wcomp_spec_vec_t,wcomp_spec_vec,wcomp_spec_t);
