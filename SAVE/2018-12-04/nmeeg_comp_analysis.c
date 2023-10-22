#define PROG_NAME "nmeeg_comp_analysis"
#define PROG_DESC "Performs principal component analysis on an EEG dataset."
#define PROG_VERS "2013-06-06"

#define nmeeg_comp_analysis_C_COPYRIGHT \
  "Copyright © 2013 by the State University of Campinas (UNICAMP)"
/* Last edited on 2023-10-21 21:57:22 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    { -pattern {PATNAME} {PATFILE} }.. \\\n" \
  "    [ -normalize ] \\\n" \
  "    [ -writeComp {COMPNAME} {COMPFILE} = {PATNAME[0]} {PATNAME[1]} .. {PATNAME[M-1]} ].. \\\n" \
  "    [ -norm {NORMNAME} = {CHNAME[0]} {CHNAME[1]} .. {CHNAME[M-1]} ].. \\\n" \
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
  " frame {V} as a vector in Cartesian space of dimension {NE}; and the" \
  " first {NE} values in each pattern frame {PAT[k]} as a similar vector {P[k]}.  It" \
  " will decompose each input vector {VI} into two components: a linear combination {VP} of" \
  " the pattern vectors\n" \
  "\n" \
  "   {VP = C[0]*P[0] + C[1]*P[1] + ··· + C[NP-1]*P[NP-1]}\n" \
  "\n" \
  " where {C[0..NP-1]} are {NP} real coefficients; and a residual" \
  " vector {VR} with {NE} elements, that is orthogonal" \
  " to all vectors {P[0..NP-1]}.   The program will then concatenate the vectors {VR} and {C}" \
  " into an output vector {VO} with {NE+NP} channels.\n" \
  "\n" \
  "  Optionally, the program may append to {VO} some additional channels, each" \
  " of them the Euclidean norm of some subset set of other channels.\n" \
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
  "   {VS = C[K[0]]*P[K[0]] + C[K[1]]*P[K[1]] + ··· + C[K[M-1]]*P[K[M-1]]}\n" \
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
  "  -pattern {PATNAME} {PATFILE}\n" \
  "    Each occurrence of this optional argument specifies the name {PATNAME} of a" \
  " pattern and the name of the file {PATFILE} that contains it.  The {PATNAME} must start with" \
  " an ASCII letter and may contain only ASCII letters, digits" \
  " 0-9, periods '.' and underscores '_'.  The {PATNAME}s must be distinct among themselves" \
  " and from other channel names.\n" \
  "\n" \
  "  -normalize\n" \
  "    This optional argument specifies that the electrode values of input patterns must be scaled after" \
  " reading so that the Euclidean norm of the electrode values is 1.  If this parameter is omitted, the" \
  " pattern frames are used as read.\n" \
  "\n" \
  "  -writeComp {COMPNAME} {COMPFILE} = {PATNAME[0]} {PATNAME[1]} .. {PATNAME[M-1]}\n" \
  "    Each occurrence of this optional argument requests the output to a" \
  " file called {COMPFILE} of an EEG dataset" \
  " consisting of a linear combination of patterns {PATNAME[0..M-1]}, with the respective" \
  " coefficients as fitted in each frame.  Each {PATNAME[K]} must be the name of a" \
  " pattern defined in a preceding \"-pattern\" argument.   The {COMPNAME} is a" \
  " string that will be written as the \"component\" field of the output" \
  " file's header; it must start with" \
  " an ASCII letter and may contain only ASCII letters, digits" \
  " 0-9, periods '.' and underscores '_'.\n" \
  "\n" \
  "  -norm {NORMNAME} = {CHNAME[0]} {CHNAME[1]} .. {CHNAME[M-1]}\n" \
  "    Each  occurrence of this optional argument specifies the name {NORMNAME} of a" \
  " new channel that is to be computed as the Euclidean norm the {M} electrode" \
  " channels {CHNAME[0]}, {CHNAME[1]}, ..., {CHNAME[M-1]}.  The latter must be" \
  " either original input channel names, or new channels defined in" \
  " previous \"-pattern\" o \"-norm\" arguments.  The norm is computed after the" \
  " patterns have been fitted to and subtracted from the original" \
  " channels.  The {NORMNAME} must start with" \
  " an ASCII letter and may contain only ASCII letters, digits" \
  " 0-9, periods '.' and underscores '_'.  The {NORMNAME}s must be distinct among themselves" \
  " and from other channel names. \n" \
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

#include <neuromat_eeg.h>
#include <neuromat_eeg_io.h>
#include <neuromat_eeg_header.h>
#include <neuromat_eeg_stats.h>

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

typedef struct norm_spec_t
  { char *normname;       /* Name of new channel defined by norm of other channels. */
    string_vec_t chname;  /* Names of channels that enter in the norm. */
  } norm_spec_t;
  
vec_typedef(norm_spec_vec_t,norm_spec_vec,norm_spec_t);

typedef struct nca_options_t
  { pat_spec_vec_t pattern;     /* Channels defined by patterns. */
    wcomp_spec_vec_t writeComp; /* Components to write separately. */
    norm_spec_vec_t norm;       /* Channels defined by norm. */
    string_vec_t delete;        /* Channels to be deleted. */
    bool_t normalize;           /* TRUE to normalize the patterns. */
  } nca_options_t;
  /* Arguments from command line. */
  
void nca_read_dataset(FILE *rd, neuromat_eeg_header_t **hP, double ***valP);
  /* Reads an EEG dataset from file {rd}, returns its sample array in {*valP} and the 
    header record {*hP}.  
    
    The sample array has the form {val[0..h.nt-1][0..h.nc-1]} where {val[t][i]}
    is the sample in channel {i} in the frame with index {t}. */
  
void nca_read_pattern(char *fname, int ne, char *chname[], bool_t normalize, double pat[]);
  /* Reads a single-frame EEG dataset from file {fname}, returns the electrode
    values in {pat[0..ne-1]}.  The dataset must have at least {ne} channels
    whose names must match {chname[0..ne-1]}. Other channels are ignored. */
        
void nca_compute_fitting_matrix(int np, int ne, double *P, double *Q);
  /* Stores into {Q} the inverse of {P*P'}.  Also checks whether the inverse is 
    valid. */

void nca_print_channel_stats(FILE *wr, int nt, int nc, int ne, double **val, char *chname[]);
  /* Gathers and prints the per-channel statistics of the data set {val[0..nt-1][0..nc-1]}.
    Assumes that the first {ne} channels are electrode voltages (real or synthetic). */

void nca_extract_patterns(int nt, int ne, double **vin, int np, double *P, double *Q, double **vot);
  /* Converts {ne} electrode potential channels from the dataset {vin}
    to {ne+np} channels in the new dataset {vot}, using the pattern matrix {P}.

    Assumes that {vin[it]} and {vot[it]} are corresponding input and
    output data frames for each {it} in {0..nt-1}; which have at least
    {ne} and {ne+np} electrode values, respectively (possibly followed
    by other channels).
    
    Also assumes that {P} is a linearized array with {np} rows and
    {ne} columns, such that {P[ip*ne + ie]} is the value of
    electrode {ie} in pattern {ip}.  Assumes that {Q} is a linearized
    matrix with {np} rows and {np} columns, the inverse of {P*P'}
    where {P'} is the transpose of {P}. 
    
    For each input frame {vin[it]}, the program considers its first {ne}
    samples to be a row vector {vi[0..ne-1]}, and computes {np}
    coefficients {vo[0..np-1]} by the product {Q*P*vi'}. These
    coefficients are stored into the output frame {vot[it]} as channels
    {ne..ne+np-1}.
    
    Then the program subtracts from the vector {vi} the fitted components, namely {P*vo'},
    and stores that residual into {vot[it]} as channels {0..ne-1}. */

void nca_recombine_patterns(int nt, int ne, int np, double **vsp, int mp, int ip[], double *P, double **vrc);
  /* Separates a specified pattern from an EEG dataset.
  
    Assumes that {vsp} is an EEG dataset with {nt} frames {vsp[0..nt-1]},
    each frame containing {ne} residual electrode readings followed by {np} pattern
    coefficients (and possibly other channels), as extracted by
    {nca_extract_patterns}. Assumes that {P} is a linearized array with
    {np} rows and {ne} columns, where {P[ip*ne + ie]} is the value of
    electrode {ie} in pattern {ip}.  Assumes that {ip[0..mp-1]} are 
    integers in the range {0..np-1}, Finally assumes that {vrc} is an 
    array of {nt} data frames, each having {nc} channels.
    
    The procedure stores into channels {0..ne-1} of each frame
    {vrc[it]} the EEG electrode potentials corresponding to the
    linear combination of rows {ip[0..mp-1]} of {P}, with their
    respective coefficients as found by {nca_extract_patterns};
    namely, sets
    
      {vrc[it][ie] = SUM{ vsp[it][ne+ip[k]]*P[ip[k]*ne + ie] : k \in 0..mp-1}}
      
    for {ie}in {0..ne-1}.  The remaining {nc-ne} channels of each frame
    {vrc[it]} are left unchanged. */
          
void nca_append_norm_channel(int nt, int nc, double **val, int mc, int ic[]);
  /* Assumes that {val} is an array of {nt} frames, each with at least {nc+1} channels
    of which the first {nc} are defined.  Also assumes that {ic[0..mc-1]} is
    a list of indices in {0..nc-1}.  Sets channel {nc} of every frame {v=val[it]}
    to the Euclidean norm of the values {v[ic[0..mc-1]}. */

void nca_delete_channels(int nt, int *ncP, int *neP, double **val, char *chname[], bool_t keep[]);
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

void nca_write_eeg_dataset(FILE *wr, int nt, int nc, double **val, neuromat_eeg_header_t *h);
  /* Writes the signals {val[0..nt-1][0..nc-1]}, comprising {nc}
    channels sampled at {nt} times, to the file "{prefix}_{tag}.txt".
    Writes the header {h} and then calls {neuromat_eeg_data_write} to
    write the data. */

nca_options_t *nca_parse_options(int argc, char **argv);
  /* Parses the command line arguments. */

string_vec_t nca_parse_chnames(argparser_t *pp);
  /* Parses a list of channel names from the command line, stopping at the next 
     keyword (any argument that starts with "-"). */

int main(int argc, char **argv);
  /* Main prog. */

int main(int argc, char **argv)
  {
    nca_options_t *o = nca_parse_options(argc, argv);
    
    /* Read the dataset: */
    neuromat_eeg_header_t *h = NULL; /* Parameters from the input file header. */
    double **vin = NULL;
    nca_read_dataset(stdin, &h, &vin);
    int nt = h->nt;
    int nc_in = h->nc;
    int ne_in = h->ne;
    
    fprintf(stderr, "input channel statistics\n");
    nca_print_channel_stats(stderr, nt, nc_in, ne_in, vin, h->chname);
    
    int np = o->pattern.ne; /* Number of pattern frames: */
    int nr = o->norm.ne;    /* Number of norm channels. */
    int nd = o->delete.ne;  /* Number of channels to be deleted. */
    
    int it, ie;

    /* Allocate the output dataset will all channels (before channel deletion): */
    int nc_tp = nc_in + np + nr; /* Number of temporary channels. */
    double **vot = notnull(malloc(nt*sizeof(double*)), "no mem"); /* Index {[0..nt-1][0..nc_tp-1]}. */
    for (it = 0; it < nt; it++) { vot[it] = notnull(malloc(nc_tp*sizeof(double)), "no mem"); }
        
    /* Allocate the list of output channel names (with space also for channels to be deleted): */
    char *chname_ot[nc_tp];
    
    /* Initialize the list of output channels with the input electrode names: */
    int nc_ot = ne_in; /* Number of channels in output file. */
    int ne_ot = ne_in; /* Number of electrodes (real or synthetic) in output file. */
    for (ie = 0; ie < ne_in; ie++) { chname_ot[ie] = txtcat(h->chname[ie], ""); }
    
    if (np > 0)
      { /* Read the patterns into {P[0..np*ne-1]}, save their names in {chname_ot[ne_in..ne_in+np-1]}: */
        demand(np <= ne_in, "too many patterns");
        double *P = rmxn_alloc(np,ne_in); /* Each row is a pattern frame. */
        char **patnames = &(chname_ot[nc_ot]); /* Section of {chname_ot} with the pattern names. */
        int ip;
        for (ip = 0; ip < np; ip++)
          { pat_spec_t *Sip = &(o->pattern.e[ip]);
            double *Pip = &(P[ip*ne_in]);
            fprintf(stderr, "reading pattern %s from file %s\n", Sip->patname, Sip->fname);
            nca_read_pattern(Sip->fname, ne_in, h->chname, o->normalize, Pip);
            /* Check for repeated names: */
            int jp = neuromat_eeg_find_channel_by_name(Sip->patname, 0, nc_ot-1, chname_ot, FALSE);
            demand(jp == -1, "repeated pattern or electrode name");
            assert(nc_ot == ne_in + ip); /* Paranoia. */
            patnames[ip] = txtcat(Sip->patname, "");
            nc_ot++;
            ne_ot++;
          }
          
        fprintf(stderr, "== debug 0 ==\n");
        int ic; for (ic = 0; ic < nc_ot; ic++) { fprintf(stderr, "channel %3d = %-4s\n", ic, chname_ot[ic]); }
          
        /* Compute the array {Q[0..np*np-1]}, the inverse of {P*P'}: */
        double *Q = rmxn_alloc(np,np); /* The array {(P*P')^-1}. */
        nca_compute_fitting_matrix(np,ne_in,P,Q);

        fprintf(stderr, "== debug 0.4 ==\n");
        for (ic = 0; ic < nc_ot; ic++) { fprintf(stderr, "channel %3d = %-4s\n", ic, chname_ot[ic]); }
        
        /* Extract the patterns: */
        nca_extract_patterns(nt, ne_in, vin, np, P, Q, vot);

        fprintf(stderr, "== debug 0.5 ==\n");
        for (ic = 0; ic < nc_ot; ic++) { fprintf(stderr, "channel %3d = %-4s\n", ic, chname_ot[ic]); }

        /* Write the fitted pattern combinations: */
        int nw = o->writeComp.ne;
        if (nw > 0)
          { int iw;
            for (iw = 0; iw < nw; iw++)
              { wcomp_spec_t *Siw = &(o->writeComp.e[iw]);  /* Spec of one output pattern combination file. */
                fprintf(stderr, "writing component %s to file %s\n", Siw->compname, Siw->fname);
                int mp = Siw->patnames.ne; /* Number of patterns to combine into this file. */
                int ip[mp]; /* Indices of patterns to combine. */
                int k;
                for (k = 0; k < mp; k++) 
                  { ip[k] = neuromat_eeg_find_channel_by_name(Siw->patnames.e[k], 0, np-1, patnames, TRUE); }
                nca_recombine_patterns(nt, ne_in, np, vot, mp, ip, P, vin);
                h->component = Siw->compname;
                FILE *wr = open_write(Siw->fname, TRUE);
                nca_write_eeg_dataset(wr, nt, nc_in, vin, h);
                fclose(wr);
              }
          }
      }

    fprintf(stderr, "== debug 1 ==\n");
    int ic; for (ic = 0; ic < nc_ot; ic++) { fprintf(stderr, "channel %3d = %-4s\n", ic, chname_ot[ic]); }

    if (nr > 0)
      { /* Add the norm channels to {vot}: */
        int ir;
        for (ir = 0; ir < nr; ir++)
          { norm_spec_t *Sir = &(o->norm.e[ir]);;
            int mc = Sir->chname.ne; /* Number of channels to be combined into this norm */
            int ic[mc]; /* Indices of channels to combine. */
            int k;
            for (k = 0; k < mc; k++) 
              { ic[k] = neuromat_eeg_find_channel_by_name(Sir->chname.e[k], 0, nc_ot-1, chname_ot, TRUE); }
            nca_append_norm_channel(nt, nc_ot, vot, mc, ic);
            int jc = neuromat_eeg_find_channel_by_name(Sir->normname, 0, nc_ot-1, chname_ot, FALSE);
            demand(jc < 0, "repeated norm, electrode, or pattern name");
            chname_ot[nc_ot] = txtcat(Sir->normname, "");
            nc_ot++;
            ne_ot++;
          }
      }

    assert(ne_ot == nc_ot);

    fprintf(stderr, "== debug 2 ==\n");
    for (ic = 0; ic < nc_ot; ic++) { fprintf(stderr, "channel %3d = %-4s\n", ic, chname_ot[ic]); }

    if (nc_in > ne_in)
      { /* Append marker/trigger (non-electrode) channels: */
        int it, ic, jc;
        for (it = 0; it < nt; it++) 
          { double *vi = vin[it];
            double *vo = vot[it];
            for (ic = ne_in, jc = nc_ot; ic < nc_in; ic++, jc++) { vo[jc] = vi[ic]; }
          }
        for (ic = ne_in; ic < nc_in; ic++)
          { char *mkname = h->chname[ic];
            int kc = neuromat_eeg_find_channel_by_name(mkname, 0, nc_ot-1, chname_ot, FALSE);
            demand(kc < 0, "collision with marker channel name");
            chname_ot[nc_ot] = txtcat(mkname, "");
            nc_ot++;
          }
      }

    assert(nc_ot == nc_tp);

    fprintf(stderr, "== debug 3 ==\n");
    for (ic = 0; ic < nc_ot; ic++) { fprintf(stderr, "channel %3d = %-4s\n", ic, chname_ot[ic]); }

    if (nd > 0)
      { /* Delete channels: */
        bool_t keep[nc_ot]; /* TRUE for channels to keep. */
        int ic, id;
        for (ic = 0; ic < nc_in; ic++) { keep[ic] = TRUE; }
        for (id = 0; id < nd; id++)
          { char *delname = o->delete.e[id];
            int ick = neuromat_eeg_find_channel_by_name(delname, 0, nc_ot-1, chname_ot, TRUE);
            demand(keep[ick], "duplicate channel in delete list");
            keep[ick] = FALSE;
          }
        nca_delete_channels(nt, &nc_ot, &ne_ot, vot, chname_ot, keep);
      }

    fprintf(stderr, "== debug 4 ==\n");
    for (ic = 0; ic < nc_ot; ic++) { fprintf(stderr, "channel %3d = %-4s\n", ic, chname_ot[ic]); }

    fprintf(stderr, "output channel statistics\n");
    nca_print_channel_stats(stderr, nt, nc_ot, ne_ot, vot, chname_ot);
    
    /* Write the output EEG dataset: */
    h->nc = nc_ot; /* Principal component strengths and triggers. */
    h->chname = chname_ot;
    h->ne = ne_ot;
    h->component = NULL;
    nca_write_eeg_dataset(stdout, nt, nc_ot, vot, h);

    /* Dont's bother to free storage: */
    return 0;
  }

void nca_delete_channels(int nt, int *ncP, int *neP, double **val, char *chname[], bool_t keep[])
  {
    int nc_in = (*ncP);
    int ne_in = (*neP);
    int ic, nc_ot, ne_ot;
    
    /* Delete the channels in all frames: */
    int it;
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

void nca_append_norm_channel(int nt, int nc, double **val, int mc, int ic[]) 
  {
    int it, k;
    for (k = 0; k < mc; k++) { demand((ic[k] >= 0) && (ic[k] < nc), "invalid channel index"); }
    for (it = 0; it < nt; it++)
      { double *v = val[it];
        assert(v != NULL);
        int k;
        double sum2 = 0;
        for (k = 0; k < mc; k++)
          { double vk = v[ic[k]];
            sum2 += vk*vk;
          }
        v[nc] = sqrt(sum2);
      }
  }
           
void nca_compute_fitting_matrix(int np, int ne, double *P, double *Q)
  {
    double *R = rmxn_alloc(np,np); /* The array {P*P'}. */
    rmxn_mul_tr(np, np, ne, P, P, R);
    rmxn_inv(np,R,Q);
    /* Check the inverse: */
    double *S = rmxn_alloc(np,np); /* Should be the identity. */
    rmxn_mul(np,np,np,R,Q,S);
    int i, j;
    for (i = 0; i < np; i++)
      { for (j = 0; j < np; j++)
          { double Sij = S[i*np + j];
            double Iij = (i == j ? 1.0 : 0.0);
            if (fabs(Sij - Iij) > 1.0e-6)
              { fprintf(stderr, "  S[%d,%d] = %24.16e\n", i, j, Sij);
                demand(FALSE, "inversion failed, aborted");
              }
          }
      }
    free(R);
    free(S);
  }
    
void nca_read_dataset(FILE *rd, neuromat_eeg_header_t **hP, double ***valP)
  {
    int nt = 0; /* Number of data frames read from all datasets. */
    
    /* Read the input file's header, define {nc,ne}: */
    int nl = 0;
    neuromat_eeg_header_t *h = neuromat_eeg_header_read(rd, 20, 600, &nl);
    int nc = h->nc; /* Number of channels in input file.*/
    demand(nc > 0, "invalid {nc} in file");
    int ne = h->ne; /* Number of electrode signals in input file. */
    demand((0 < ne) && (ne <= nc), "invalid {ne} in file");
    demand (h->chname != NULL, "missing channel names in input file header");
    /* int ic; for (ic = 0; ic < nc; ic++) { fprintf(stderr, "channel %3d = %-4s\n", ic, h->chname[ic]); } */
        
    /* Read the EEG data, define {nt} (ignore header {nt} if any): */
    double **val = neuromat_eeg_data_read(rd, 0, h->nt, nc, &nl, &nt);

    fprintf(stderr, "read %d data frames\n", nt);
    fprintf(stderr, "input file has %d channels including %d electrodes\n", nc, ne);
    demand(nt > 0, "no frames to process");
    
    (*hP) = h;
    (*valP) = val;
  }
    
void nca_read_pattern(char *fname, int ne, char *chname[], bool_t normalize, double pat[])
  {
    FILE *rd = open_read(fname, TRUE); 
    neuromat_eeg_header_t *h = NULL;
    double **val = NULL;
    nca_read_dataset(rd, &h, &val);
    demand(h->ne == ne, "electrodes counts in pattern and file do not match");
    demand(h->nt == 1, "pattern should have a single data frame");
    demand (h->chname != NULL, "missing channel names in pattern header");
    double sum2 = 0;
    int ie;
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
      }
    free(val);
    free(h);
    fclose(rd);
  }
    
void nca_print_channel_stats(FILE *wr, int nt, int nc, int ne, double **val, char *chname[])
  {
    /* Compute the channel means and variances over the selected frames: */
    double vavg[nc], vvar[nc], vmin[nc], vmax[nc];
    double vmin_min, vmax_max, vdev_max;  /* Maximum variance among electrode channels. */
    neuromat_eeg_stats_per_channel(nt, nc, val, FALSE, vavg, vvar, vmin, vmax, ne, &vmin_min, &vmax_max, &vdev_max);
    fprintf(wr, "--- input channel statistics ---\n");
    neuromat_eeg_stats_per_channel_print(wr, nc, chname, vavg, vvar, vmin, vmax);
    fprintf(wr, "global electrode minimum = %8.5f\n", vmin_min);
    fprintf(wr, "global electrode maximum = %8.5f\n", vmax_max);
    fprintf(wr, "maximum electrode dev = %8.5f\n", vdev_max);
    fprintf(wr, "\n");
  }

void nca_extract_patterns(int nt, int ne, double **vin, int np, double *P, double *Q, double **vot)
  {
    demand(ne >= 0, "invalid {ne}");
    demand(np >= 0, "invalid {np}");
    demand(np <= ne, "more components than electrodes");
    int it;
    for (it = 0; it < nt; it++)
      { double *vi = vin[it]; /* Original electrode values {vi[0..ne-1]}. */
        double *cf = &(vot[it][ne]); /* Coefficients of each pattern {cf[0..np-1]}. */
        double *vo = vot[it]; /* Residual is {vo[0..ne-1]}. */
        /* Map channels {valt[0..ne-1]} by first {np} rows of {P} to {vpct[0..np-1]}: */
        double b[np];
        rmxn_map_col(np, ne, P, vi, b);
        rmxn_map_col(np, np, Q, b, cf);
        /* Compute residual: */
        double vr[ne]; /* Reconstructed frame from fitted patterns. */
        rmxn_tr_mul(np, ne, 1, P, cf, vr);
        rn_sub(ne, vi, vr, vo);
      }
  }

void nca_recombine_patterns(int nt, int ne, int np, double **vsp, int mp, int ip[], double *P, double **vrc)
  {
    demand((0 <= np) && (np <= ne), "invalid {np,ne}");
    int it;
    for (it = 0; it < nt; it++)
      { double *cf = &(vsp[it][ne]); /* Fitted pattern coefficients. */
        double *vr = vrc[it];        /* Reconstructed electrodes. */
        /* Combine the requested components: */
        rn_zero(ne, vr);
        int k;
        for (k = 0; k < mp; k++)
          { int ipk = ip[k];
            demand((0 <= ipk) && (ipk < np), "invalid pattern index");
            double *Pi = &(P[ipk*ne]); /* Selected pattern. */
            double cfi = cf[ipk]; /* Its coefficient. */
            rn_mix_in (ne, cfi, Pi, vr);
          }
      }
  }

void nca_write_eeg_dataset(FILE *wr, int nt, int nc, double **val, neuromat_eeg_header_t *h)
  {
    assert(h->nt == nt);
    assert(h->nc == nc);
    neuromat_eeg_header_write(wr, h);
    neuromat_eeg_data_write(wr, nt, nc, val, 0, nt-1, 1);
    fflush(wr);
  }

nca_options_t *nca_parse_options(int argc, char **argv)
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
      int np = 0;
      while (argparser_keyword_present(pp, "-pattern"))
        { char *patname = argparser_get_next_non_keyword(pp);
          char *fname = argparser_get_next_non_keyword(pp);
          pat_spec_vec_expand(&(o->pattern), np);
          o->pattern.e[np] = (pat_spec_t){ .patname = patname, .fname = fname };
          np++;
        }
      pat_spec_vec_trim(&(o->pattern), np);
    }
    
    { o->norm = norm_spec_vec_new(10);
      int nr = 0;
      while (argparser_keyword_present(pp, "-norm"))
        { char *normname = argparser_get_next_non_keyword(pp);
          argparser_get_keyword_next(pp, "=");
          string_vec_t chname = nca_parse_chnames(pp);
          norm_spec_vec_expand(&(o->norm), nr);
          o->norm.e[nr] = (norm_spec_t){ .normname = normname, .chname = chname };
          nr++;
        }
      norm_spec_vec_trim(&(o->norm), nr);
    }

    { o->writeComp = wcomp_spec_vec_new(10);
      int nw = 0;
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

    { o->delete = string_vec_new(10);
      int nd = 0;
      while (argparser_keyword_present(pp, "-delete"))
        { string_vec_t chname = nca_parse_chnames(pp);
          int id;
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
    int nc = 0;
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
vec_typeimpl(norm_spec_vec_t,norm_spec_vec,norm_spec_t);
vec_typeimpl(wcomp_spec_vec_t,wcomp_spec_vec,wcomp_spec_t);
