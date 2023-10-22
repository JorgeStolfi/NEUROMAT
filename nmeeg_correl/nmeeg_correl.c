#define PROG_NAME "nmeeg_correl"
#define PROG_DESC "Performs principal component analysis on an EEG dataset."
#define PROG_VERS "2013-06-06"

#define nmeeg_correl_C_COPYRIGHT \
  "Copyright © 2013 by the State University of Campinas (UNICAMP)"
/* Last edited on 2023-10-22 09:42:49 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -capType {CAP_TYPE} \\\n" \
  "    [ -select {SNAME} {SVAL} ] \\\n" \
  "    [ -normAvg {NORM_AVG} ] \\\n" \
  "    [ -normVar {NORM_VAR} ] \\\n" \
  "    [ -minMag {MINMAG} ] \\\n" \
  "    [ -maxComps {MAXCOMPS} ] \\\n" \
  "    -outPrefix {OUT_PREFIX} \\\n" \
  "    " argparser_help_info_HELP " \\\n" \
  "    < {INFILE}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program reads from standard input one or more EEG datasets" \
  " and determines the principal components of their electrode signals.\n" \
  "\n" \
  "  The program assumes that every data line (/data frame/) in the input file contains" \
  " simultaneous samples from some fixed number {NC} of data channels, of which" \
  " the first {NE} are electrode potentials, while the rest are phase indicators (/markers/)" \
  " or other unspecified signals.\n" \
  "\n" \
  "  The program computes the covariance matrix {M} of the {NE} electrode" \
  " signals, after optionally subtracting from each signal its average value" \
  " and normalizing its variance to 1, and then" \
  " extracts the {NV} principal components (PC) from it.  Each PC is an" \
  " eigenvector of {M}, with unit Euclidean norm.  The corresponding" \
  " eigenvalue is the square of the RMS amplitude of that component.  The" \
  " number {NV} is at most {NE}; it is less if the numerical eigendecomposition" \
  " algorithm fails.\n" \
  "\n" \
  "  The program can be instructed to use only the" \
  " data frames where a particular marker channel has a specific value.\n" \
  "\n" \
  "INPUT\n" \
  "  The input file should contain one or more EEG datasets.  Each dataset" \
  " should have one or more header lines that specify (among other parameters) the" \
  " number of channels {NC}, the number of electrodes {NE}, and the number of data" \
  " frames {NT}.  All datasets must have the same {NE} and {NC} parameters.  After" \
  " the header lines there must be {NT} data lines, each containing {NC} floating-point" \
  " numbers.\n" \
  "\n" \
  "OUTPUTS\n" \
  "  The program writes to a file called \"{OUT_PREFIX}dst_cov.txt\" the" \
  " distance and covariance information for all pairs of electrodes.  Each" \
  " line has fields \"{I} {J} {D2IJ} {CVIJ} {CVII} {CVJJ}\" where {D2IJ} is" \
  " the square of the distance between the 3D positions of electrodes {I} and" \
  " {J}, {CVIJ} is the covariance of those elecrodes, and {CVII} and {CVJJ} are" \
  " their variances (after channel normalization).\n" \
  "\n" \
  "  Each principal component (eigenvector) is written out as a single-frame dataset, with" \
  " the original {NC} channels and {NE} electrodes, to a file" \
  " called \"{OUT_PREFIX}P{IV}_eig.txt\" where {IV} is the component's" \
  " index in {0..NV-1}.  The channels are scaled so that their RMS value is 1.\n" \
  "\n" \
  "  Output files will have some header records before the first data frame.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -capType {CAP_TYPE}\n" \
  "    This mandatory argument specifies the cap type, which determines" \
  " the count, names, and positions of the electrodes.  The possible" \
  " values of {CAP_TYPE} are defined" \
  " by {neuromat_eeg_get_channel_names} in {neuromat_eeg.h}, and" \
  " include \"R20\", \"R128\", \"R129\", and \"FN3\".\n" \
  "\n" \
  "  -select {SNAME} {SVAL}\n" \
  "    This optional argument specifies the name {SNAME} of a marker channel that" \
  " may be used to select which frames to use when computing the covariance" \
  " matrix {M}; namely, only those frames where that channel has" \
  " the specified value {SVAL}.  The marker channel's position in each data frame is obtained" \
  " from the list of channel names found in the input file's header.  If this" \
  " argument is omitted, all input frames are used.\n" \
  "\n" \
  "  -normAvg {NORM_AVG}\n" \
  "    This optional argument asks for average value normalization.  If {NORM_AVG} is \"T\" or 1, the" \
  " program computes the mean value of each electrode signal and" \
  " subtracts it before computing the covariances.  If {NORM_AVG} is \"F\" or 0, the" \
  " program assumes that the mean of each signal is zero.  If omitted, the program" \
  " assumes \"-normAvg T\".\n" \
  "\n" \
  "  -normVar {NORM_VAR}\n" \
  "    This optional argument asks for variance normalization.  If {NORM_VAR} is \"T\" or 1, the" \
  " program divides each electrode signal" \
  " by its standard deviation before computing the covariances.  If {NORM_VAR} is \"F\" or 0, the" \
  " program does not do any scaling.  If omitted, the program assumes \"-normVar F\".  Note" \
  " that if {NORM_AVG} is false the deviation is computed assuming that the mean is 0.\n" \
  "\n" \
  "  -minMag {MINMAG}\n" \
  "    This optional argument specifies the minimum RMS strength of a principal" \
  " component (euclidean norm of the electrode sample vector) in microvolts (uV).  If" \
  " omitted, the program assumes {MINMAG=1.0}.\n" \
  "\n" \
  "  -maxComps {MAXCOMPS}\n" \
  "    This optional argument specifies the maximum number of principal" \
  " components to output.  If it is negative or omitted, the program outputs" \
  " all principal coponents that satisy the \"-minMag\" criterion.\n" \
  "\n" \
  "  -outPrefix {OUT_PREFIX}\n" \
  "    This mandatory argument specifies a common prefix for output" \
  " file names (see the DESCRIPTION section above).  If the {OUT_PREFIX} string" \
  " includes any slashes, the corresponding directories must exist" \
  " and be writable.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  neuromat_eeg_plot_signals.sh(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2013-06-04 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2013-06-06 J. Stolfi: Converted to {argparser.h}.\n" \
  "  2013-06-09 J. Stolfi: Added anomaly channel, made \"-writeComps\" optional.\n" \
  "  2013-11-24 J. Stolfi: allowed multiple input datasets, split off the decomposition.\n" \
  "  2013-11-29 J. Stolfi: added \"-maxComps\" option, removed vestiges of anomaly channel.\n" \
  "  2013-11-30 J. Stolfi: added \"-zeroMean\" option.\n" \
  "  2021-08-23 J. Stolfi: replaced \"-zeroMean\" by \"-normAvg\", added \"-normVar\".\n" \
  "  2023-10-22 J. Stolfi: added \"-capType\" parameter.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " nmeeg_correl_C_COPYRIGHT ".\n" \
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

#include <rmxn.h>
#include <rn.h>
#include <r2.h>
#include <r3.h>
#include <fget.h>
#include <argparser.h>
#include <sym_eigen.h>
#include <affirm.h>
#include <jsstring.h>
#include <jsfile.h>
#include <jsmath.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_io.h>
#include <neuromat_eeg_geom.h>
#include <neuromat_eeg_header.h>
#include <neuromat_eeg_channel_stats.h>
#include <neuromat_eeg_pca.h>

typedef struct nec_options_t
  {
    char *capType;        /* Type of EEG cap, like "R20", "R128" etc. */
    char *select_sname;   /* Name of marker channel, or {NULL}. */
    double select_sval;   /* Marker value indicateing valid data, or {NAN}. */
    bool_t normAvg;       /* TRUE to subtract the mean from each electrode signal. */
    bool_t normVar;       /* TRUE to scale each signal by its standard deviation from 0. */
    double minMag;        /* Minimum strength of a principal component. */
    int maxComps;         /* Max number of components to output, or negative for "all". */
    char *outPrefix;      /* Prefix for all output filenames. */
  } nec_options_t;
  /* Arguments from command line. */
  
double **nec_read_file(FILE *rd, nec_options_t *o, neuromat_eeg_header_t **hP);
  /* Reads one or more EEG datasets from file {rd}, returns the concatenation 
    of all sample arrays, and also builds a combined header record {*hP} for them.  
    
    The result is an array {val[0..h.nt-1][0..h.nc-1]} where {val[t][i]}
    is the sample in channel {i} in the frame with index {t}. This
    dataset may contain only a subset of the data, selected by the
    parameters {o->select_sname,o->select_sval}. */
  
neuromat_eeg_header_t *nec_read_header(FILE *rd, int *nlP);
  /* Reads a single dataset from file {rd}, returns its header record. 
    Returns NULL if there is no header. */
  
vec_typedef(nec_frame_vec_t,nec_frame_vec,double*);
  /* A vector of pointers to frames. */
  
void nec_read_data_frames(FILE *rd, int nt_read, int nc, int *nlP, int *nfP, int ict, double valt, int *ntP, nec_frame_vec_t *fvP);
  /* Reads from file {rd} a set of {nt_read} data frames, appends some of them to {*fvP}.  
    
    If {ict} is non-negative, the procedure considers only data frames
    where channel {ict} has value {valt}. Otherwise all frames read
    from {rd} are selecetd.
    
    The selected frames are stored in {fvP.e[nt..nt+st-1][0..nc-1]},
    where {nt} is the input value of {*ntP} and {st} is the number of
    frames that were read and selected. The procedure updates {*ntP} to
    {nt+sp}.
    
    The procedure assumes that {*nlP} is the number of lines already
    read from {rd} (including header and comment lines), and increments
    {*nlP} with the number of lines actually read (including any
    '#'-comments and header lines).
    
    The procedure also assumes that {*nfP} is the number of data frames
    already read from {rd} (excluding header and comment lines), and
    increments {*nfP} by the number of data frames that it reads (0 or 1). */

void nec_channel_stats(FILE *wr, neuromat_eeg_header_t *h, double **val, bool_t normAvg, bool_t normVar, double vshift[]);
  /* Gathers and prints the per-channel statistics of the data set {val[0..h.nt-1][0..h.nc-1]}.
    Returns in {vshift[0..h.nc-1]} the average of each channel.   If {normAvg} is faslse,
    computes the covariances assuming as if the expected value of each electrode signal is zero. */

void nec_compute_principal_components
  ( int ne, 
    double *Cv,  
    double minMag, 
    int *nvP, 
    double Ev[], 
    double emag[]
  );
  /* Finds a set of principal components for a set of {ne}
    signals, given its linearized covariance matrix 
    {Cv[0..ne*ne-1]}.

    The principal components will be the eigenvalues of the covariance
    matrix whose RMS strength in the dataset (the square root of the
    corresponding eigenvalue) is at least {minMag}. Returns in {*nvP}
    the number {nv} of components that were chosen. The components are
    stored into the first {nv} rows of the linearized {ne} by {ne}
    matrix {Ev[0..ne*ne-1]}. Each row has unit Euclidean norm and the
    rows are pairwise orthogonal.
    
    Also returns in {emag[i]} the RMS magnitude of the component {i} in
    the dataset. The rows of {Ev} are sorted by decreasing magnitude. */

void nec_write_principal_components
  ( char *prefix,
    int nv, 
    int ne, 
    double Ev[], 
    double emag[], 
    char **pc_name, 
    neuromat_eeg_header_t *h
  );
  /* Writes to disk the first {nv} principal components of a dataset
    with {ne} electrodes.
    
    Assumes that the eigenvectors of the principal components are stored
    in the first {nv} rows of {Ev}, which is assumed to be
    an {ne} by {ne} matrix linearized by rows. Each eigenvector
    is assumed to have unit norm; row {i} will be scaled 
    by {emag[i]} before wrtiting.
    
    Each scaled eigenvector is written as a EEG dataset file called
    "{prefix}{pc_name[i]}_eig.txt" with {ne} channels and a single
    data frame. The header of the file is derived from the input file
    header's {h}. */
    
void nec_write_correlation_distance_plot(char *prefix, int ne, double *Cv, r2_t pos2D[]);
/* Assumes that {Cv[0..ne*ne-1]} is an {ne} by {ne} matrix linearized by rows,
  whose element {[i,j]} is the covariance of electrodes {i} and {j}.
  Also assumes that {pos2D[i]} is the nominal position of electrode {i}
  in the 2D schematic where the head is the unit circle. 
  
  Writes to a file called "{prefix}dst_cov.txt" a set of lines,
  one for each pair of electrodes {i,j} in {0..ne-1}.
  Each line has the indices {i} and {j}, the nominal 3D distance
  squared between them {neuromat_eeg_geom_3D_dist_sqr(pos2D[i],pos2D[j])}, 
  the covariance {Cv[i,j]} and the variances {Cv[i,i]} and {Cv[j,j]}. */

char **nec_principal_component_names(int nv, char *tag);
  /* Returns a new list {pc_name[0..nv-1]} of suitable names for 
    principal components.  They are newly created strings of the form
    "{tag}{NNN}" where {NNN} is a 3-digit component index from 0 to
    {nv-1}. */

void nec_write_eeg_eigenvector(char *prefix, char *tag, int ne, double evec[], neuromat_eeg_header_t *h);
  /* Writes the signal eigenvector {evec[0..ne-1]}, comprising {ne}
    electrode potentials, to the file "{prefix}{tag}.txt". Writes the
    header {h} then calls {neuromat_eeg_frame_write} to write the
    data. */

void nec_write_eeg_dataset(char *prefix, char *tag, int nt, int nc, double **val, neuromat_eeg_header_t *h);
  /* Writes the signals {val[0..nt-1][0..nc-1]}, comprising {nc}
    channels sampled at {nt} times, to the file "{prefix}{tag}.txt".
    Writes the header {h} and then calls {neuromat_eeg_data_write} to
    write the data. */

nec_options_t *nec_parse_options(int argc, char **argv);
  /* Parses the command line arguments. */

int main(int argc, char **argv);
  /* Main prog. */

int main(int argc, char **argv)
  {
    nec_options_t *o = nec_parse_options(argc, argv);

    fprintf(stderr, "output files are named %s*\n", o->outPrefix);
    
    /* Read the dataset and collect selected frames: */
    neuromat_eeg_header_t *h = NULL; /* Parameters from the input file header. */
    double **val = nec_read_file(stdin, o, &h);
    int nc = h->nc;
    int ne = h->ne;
    int nt = h->nt;
    
    /* Compute and print per-channel statistics: */
    double vshift[nc]; /* Channel averages. */
    nec_channel_stats(stderr, h, val, o->normAvg, o->normVar, vshift);

    /* Compute the electrode covariance matrix: */
    fprintf(stderr, "computing covariances...\n");
    double *Cv = neuromat_eeg_channel_stats_covariance_matrix(nt, ne, val, vshift, NULL);
    
    /* Get the nominal 2D electrode positions: */
    fprintf(stderr, "getting electrode positions...\n");
    r2_t *pos2D = neuromat_eeg_geom_get_schematic_2D_points_by_name(o->capType, ne, h->chname);
    
    /* Write the correlation-distance plot: */
    nec_write_correlation_distance_plot(o->outPrefix, ne, Cv, pos2D);

    /* Compute the principal components and choose their names: */
    int nv; /* Number of principal components (significant eigenvectors): */
    double *Ev = rmxn_alloc(ne,ne); /* The first {nv} rows are the selected principal components. */
    double emag[ne]; /* Square roots of eigenvalues. */
    nec_compute_principal_components(ne, Cv, o->minMag, &nv, Ev, emag);
    /* Trim excess components: */
    if ((o->maxComps >= 0) && (nv > o->maxComps)) { nv = o->maxComps; }
    
    /* Create the coponent names: */
    char **pc_name = nec_principal_component_names(nv, "P"); /* Index {0..nv-1}. */
    
    free(Cv);

    /* Write the eigenvalues as single-frame datasets: */
    nec_write_principal_components(o->outPrefix, nv, ne, Ev, emag, pc_name, h);

    /* Dont's bother to free storage: */
    return 0;
  }

double **nec_read_file(FILE *rd, nec_options_t *o, neuromat_eeg_header_t **hP)
  {
    int nd = 0; /* Number of datasets in input file. */
    int nf = 0; /* Number of data frames read from all datasets. */
    int nl = 0; /* Total number of lines read from all datasets. */
    int nt = 0; /* Numberof frames read and selected in all datasets.*/
    
    int ict = -1;  /* Index of marker channel, or -1 if no marker-selection: */
    double valt = o->select_sval; /* Value of marker channel to select.*/
    
    neuromat_eeg_header_t *h = neuromat_eeg_header_new(); /* Combined header. */
    
    nec_frame_vec_t fv = nec_frame_vec_new(0); /* Frames read are {fv.e[0..nt-1]}. */
   
    while (TRUE)
      { /* Read the input file's header, define {nc,ne}: */
        neuromat_eeg_header_t *hd = nec_read_header(rd, &nl);
        if (hd == NULL) { break; }
        if (nd == 0)
          { h->nc = hd->nc;
            h->ne = hd->ne;
            fprintf(stderr, "input frames have %d channels including %d electrodes\n", h->nc, h->ne);
            demand(h->nc > 0, "invalid {nc} in file header");
            demand((0 < h->ne) && (h->ne <= h->nc), "invalid {ne} in file header");
            h->fsmp = hd->fsmp;
            fprintf(stderr, "sampling frequency is %f Hz\n", h->fsmp);
            demand((! isnan(h->fsmp)) && (h->fsmp > 0), "invalid {fsmp} in file header");
            h->chname = hd->chname;
            demand (h->chname != NULL, "missing channel names in input file header");
            if (o->select_sname != NULL)
              { /* Find marker channel index {ict}, or -1: */
                ict = 0;
                while ((ict < h->nc) && (strcmp(o->select_sname, h->chname[ict]) != 0)) { ict++; }
              }
            fprintf(stderr, "selecting frames where channel %d is %.10g\n", ict, valt);
            demand(ict < h->nc, "invalid marker channel name");
          }
        else
          { demand(hd->nc == h->nc, "inconsistent channel counts");
            demand(hd->ne == h->ne, "inconsistent electrode counts");
            demand(hd->fsmp == h->fsmp, "inconsistent sampling frequencies");
            int ic;
            for (ic = 0; ic < h->nc; ic++) 
              { demand(strcmp(hd->chname[ic], h->chname[ic]) == 0, "inconsistent channel names"); }
          }
        nec_read_data_frames(rd, hd->nt, h->nc, &nl, &nf, ict, valt, &nt, &fv);
        nd++;
      }
    fprintf(stderr, "\n");
    fprintf(stderr, "read %d datasets in input file\n", nd);
    fprintf(stderr, "total %d lines %d data frames in input file\n", nl, nf);
    fprintf(stderr, "accepted total %d data frames\n", nt);
    demand(nt > 0, "no frames to process");

    h->nt = nt;    
    (*hP) = h;
    nec_frame_vec_trim(&fv, nt);
    return fv.e;
  }
    
neuromat_eeg_header_t *nec_read_header(FILE *rd, int *nlP)
  {
    fget_skip_formatting_chars(rd);
    int c = fgetc(rd);
    if (c == EOF) { return NULL; }
    ungetc(c, rd);
    neuromat_eeg_header_t *h = neuromat_eeg_header_read(rd, 20, 600, nlP);
    return h;
  }

void nec_read_data_frames(FILE *rd, int nt_read, int nc, int *nlP, int *nfP, int ict, double valt, int *ntP, nec_frame_vec_t *fvP)
  {
    int nf0 = (*nfP);
    int nl0 = (*nlP);
    int nt0 = (*ntP);
    
    int nt = nt0;
    int it;
    for (it = 0; it < nt_read; it++)
      { double *frm = notnull(malloc(nc*sizeof(double)), "no mem");
        int nr = neuromat_eeg_frame_read(rd, nc, frm, nlP, nfP);
        if (nr == 0) { break; }
        assert(nr == 1);
        if ((ict < 0) || (frm[ict] == valt))
          { /* Take this frame: */
            nec_frame_vec_expand(fvP, nt);
            fvP->e[nt] = frm; nt++;
          }
      }
    (*ntP) = nt;
  
    fprintf(stderr, "read %d lines %d data frames in dataset\n", (*nlP)-nl0, (*nfP)-nf0);
    fprintf(stderr, "accepted %d data frames in dataset\n", nt - nt0);
    return;
  }
    
void nec_channel_stats(FILE *wr, neuromat_eeg_header_t *h, double **val, bool_t normAvg, bool_t normVar, double vshift[])
  {
    int nc = h->nc;
    int ne = h->ne;
    int nt = h->nt;

    /* Compute the channel means and variances over the selected frames: */
    neuromat_eeg_channel_stats_t *st = neuromat_eeg_channel_stats_new(nc);
    neuromat_eeg_channel_stats_t *stg = neuromat_eeg_channel_stats_new(1);
    double eps = 0.01; /* Assumed uncertainty of measurement (µV). */
    neuromat_eeg_channel_stats_gather_all(nt, nc, val, NULL, eps, st, ne, stg);
    fprintf(stderr, "  --- channel statistics ---\n");
    neuromat_eeg_channel_stats_print_all(stderr, 2, nc, h->chname, FALSE, st, ne, stg);
    fprintf(stderr, "\n");

    /* Choose the shifts: */
    for (int ic = 0; ic < nc; ic++) { vshift[ic] = (normAvg ? st[ic].avg : 0.0); }

    free(st);
    free(stg);
  }
    
void nec_compute_principal_components
  ( int ne, 
    double *Cv, 
    double minMag, 
    int *nvP, 
    double Ev[], 
    double emag[]
  )
  {
    /* Find eigenvectors and their magnitudes: */
    fprintf(stderr, "computing eigenvectors...\n");
    int nv = neuromat_eeg_pca_eigen_decomp(ne, Cv, minMag, Ev, emag); /* {nv} is num of eigens actually computed. */
    assert(nv <= ne);
    if (nv < ne) { fprintf(stderr, "found only %d out of %d significant eigencomponents\n", nv, ne); }
    
    (*nvP) = nv;
  }
    
void nec_write_principal_components
  ( char *prefix,
    int nv, 
    int ne, 
    double Ev[], 
    double emag[], 
    char **pc_name, 
    neuromat_eeg_header_t *h
  )
  {
    fprintf(stderr, "writing principal components...\n");
    int kv;
    for (kv = 0; kv < nv; kv++)
      { 
        char *prefk = NULL;  /* Prefix for files of PCA component {kv}. */
        asprintf(&prefk, "%s%s_", prefix, pc_name[kv]);
        
        /* Check that eigenvector has unit norm: */
        double *Evk = &(Ev[kv*ne]);
        double Evk_norm = rn_norm(ne, Evk);
        affirm(fabs(Evk_norm - 1.0) < 1.0e-6, "eigenvector is not normalized");

        /* Scale the eigenvector by its magnitude: */
        double Evk_mag[ne]; /* Eigenvector scaled by root of eigenvalue: */
        rn_scale(ne, emag[kv], Evk, Evk_mag);

        /* Write the magnitude-scaled eigenvector as single-frame dataset: */
        h->nt = 1;
        h->nc = ne; /* Exclude trigger channels and anomaly, so {ne} not {nc}. */
        h->chname = h->chname; /* Will use only the first {ne} entries. */
        h->ne = ne; 
        h->component = pc_name[kv];
        nec_write_eeg_eigenvector(prefk, "eig", ne, Evk_mag, h);

        free(prefk);
      }
  }

void nec_write_correlation_distance_plot(char *prefix, int ne, double *Cv, r2_t pos2D[])
  {
    fprintf(stderr, "writing plot of correlation x distance...\n");
    char *fname = NULL;
    asprintf(&fname, "%sdst_cov.txt", prefix);
    FILE *wr = open_write(fname, TRUE);
    
    r3_t srad =  (r3_t){{ 150.0, 180.0, 150.0 }};
    
    int i, j;
    for (i = 0; i < ne; i++)
      { for (j = 0; j < i; j++)
          { /* Covariances of electrodes {i,j}: */
            double Cij = Cv[i*ne + j];
            double Cji = Cv[j*ne + i];
            assert(Cij == Cji);
            double Cii = Cv[i*ne + i];
            double Cjj = Cv[j*ne + j];
            
            /* Distance of electrodes {i,j}: */
            r3_t pi = neuromat_eeg_geom_3D_from_2D(&(pos2D[i]), NULL, &srad);
            r3_t pj = neuromat_eeg_geom_3D_from_2D(&(pos2D[j]), NULL, &srad);
            double d2ij = r3_dist_sqr(&pi, &pj);
            
            fprintf(wr, "%5d %5d %11.8f %+15.8e %15.8e %15.8e\n", i, j, d2ij, Cij, Cii, Cjj);
          }
      }
    fclose(wr);
  }

void nec_write_eeg_eigenvector(char *prefix, char *tag, int ne, double evec[], neuromat_eeg_header_t *h)
  {
    assert(h->nt == 1);
    assert(h->nc == ne);
    assert(h->ne == ne);

    char *fname = NULL;
    asprintf(&fname, "%s%s.txt", prefix, tag);
    FILE *wr = open_write(fname, TRUE);
    neuromat_eeg_header_write(wr, h);

    neuromat_eeg_frame_write(wr, ne, evec);
    fclose(wr);
    free(fname);
  }

void nec_write_eeg_dataset(char *prefix, char *tag, int nt, int nc, double **val, neuromat_eeg_header_t *h)
  {
    assert(h->nt == nt);
    assert(h->nc == nc);
    
    char *fname = NULL;
    asprintf(&fname, "%s%s.txt", prefix, tag);
    FILE *wr = open_write(fname, TRUE);
    neuromat_eeg_header_write(wr, h);
    neuromat_eeg_data_write(wr, nt, nc, val, 0, nt-1, 1);
    fclose(wr);
    free(fname);
  }

char **nec_principal_component_names(int nv, char *tag)
  {
    demand(nv >= 0, "invalid {nv}");
    char **pc_name = notnull(malloc(nv*sizeof(char*)), "no mem");
    int iv;
    for (iv = 0; iv < nv; iv++)
      { char *lab = NULL;
        asprintf(&lab, "%s%03d", tag, iv);
        pc_name[iv] = lab;
      }
    return pc_name;
  }

#define nec_MAX_MARK_VAL 10000.0
  /* Max marker channel value. */

#define nec_MAX_MAXCOMPS 10000
  /* Max value of "-maxComps" parameter. */

nec_options_t *nec_parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    /* Allocate the command line argument record: */
    nec_options_t *o = notnull(malloc(sizeof(nec_options_t)), "no mem");

    /* Parse keyword parameters: */

    argparser_get_keyword(pp, "-capType");
    o->capType = argparser_get_next_non_keyword(pp);
    
    if (argparser_keyword_present(pp, "-select"))
      { o->select_sname = argparser_get_next_non_keyword(pp);
        o->select_sval = argparser_get_next_double(pp, -nec_MAX_MARK_VAL, +nec_MAX_MARK_VAL);
      }
    else
      { o->select_sname = NULL; o->select_sval = NAN; }

    if (argparser_keyword_present(pp, "-normAvg"))
      { o->normAvg = argparser_get_next_bool(pp);
      }
    else
      { o->normAvg = TRUE; }

    if (argparser_keyword_present(pp, "-normVar"))
      { o->normVar = argparser_get_next_bool(pp);
      }
    else
      { o->normVar = TRUE; }

    if (argparser_keyword_present(pp, "-minMag"))
      { o->minMag = argparser_get_next_double(pp, 0, +INF);
      }
    else
      { o->minMag = 1.0; }

    if (argparser_keyword_present(pp, "-maxComps"))
      { o->maxComps = (int)argparser_get_next_int(pp, 0, nec_MAX_MAXCOMPS);
      }
    else
      { o->maxComps = -1; }

    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next_non_keyword(pp);

    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }

vec_typeimpl(nec_frame_vec_t,nec_frame_vec,double*);
