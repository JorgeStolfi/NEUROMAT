#define PROG_NAME "nmeeg_animate"
#define PROG_DESC "Performs principal component analysis on an EEG dataset."
#define PROG_VERS "2013-06-06"

#define nmeeg_animate_C_COPYRIGHT \
  "Copyright © 2013 by the State University of Campinas (UNICAMP)"
/* Last edited on 2023-10-21 21:57:45 by stolfi */
    
/* !!! Add options to control the time bar and marker dots !!! */

/* !!! Add option to specify {vmax} !!! */

/* !!! Paint all marker channels on timebar? !!! */

/* !!! Colorize marker channels discontinuously around zero. !!! */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -basis {BASPREFIX} {BASTYPE} \\\n" \
  "    [ -skip {NSKIP} ] \\\n" \
  "    [ -read {NREAD} ] \\\n" \
  "    [ -step {FSTEP} ] \\\n" \
  "    -outPrefix {OUTPREFIX} \\\n" \
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
  "  The program reads from standard input an EEG dataset" \
  " and outputs a set of images, one for each data frame, showing the electrode intensities" \
  " as colored spots located around each electrode, overlaying a schematic diagram of" \
  " the head.  These images can be converted into an animation of the dataset.\n" \
  "\n" \
  "  The program assumes that every data line (/data frame/) in the input file contains" \
  " simultaneous samples from some fixed number {NC} of data channels, of which" \
  " the first {NE} are electrode potentials, while the rest are phase indicators (/markers/)" \
  " or other unspecified signals.\n" \
  "\n" \
  "  The program reads a /mask image/ {MSK}, a monochromatic image that defines the schematic outline" \
  " of the head; an /electrode mask/ {ELP} that shows the schematic position of each electrode; and" \
  " an /image basis/, a set of {NE} single-channel image {BIM[0..NE-1]} that define" \
  " how each electrode value is to be rendered.\n" \
  "\n" \
  "  For each input frame {VI}, the program combines the images {BIM[0..NE-1]} pixel by pixel" \
  " into a single monochromatic image {IMG}, using the electrode values as coefficients.  That is,\n" \
  "\n" \
  "   {IMG = VI[0]*BIM[0] + VI[1]*BIM[1] + ··· + VI[NE-1]*BIM[NE-1]}\n" \
  "\n" \
  " The program will then clip the image {IMG} by the mask {MSK}, convert each pixel" \
  " to an RGB color value according to a built-in false-color scale, overlay the electrode" \
  " position image {ELP} to it, and output {IMG} to a" \
  " separate file.\n" \
  "\n" \
  "FILE FORMATS\n" \
  "  The basis, mask, and electrode position files must be float-valued images" \
  " as written by {float_image_write}. All images must have a single channel and the same" \
  " size (rows and columns) in pixels.\n" \
  "\n" \
  "  The basis images may have positive and negative values.\n" \
  "\n" \
  "  The mask image samples must be 0 in the background, 1 inside the" \
  " schematic head.  The electrode position image must be 0 for" \
  " transparent, 1 for opaque.  Both should be antialiased (between" \
  " 0 and 1) along the boundary.\n" \
  "\n" \
  "  Output image files will be in PNG format.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -basis {BASPREFIX} {BASTYPE}\n" \
  "    These mandatory arguments specify the file names of the basis" \
  " images, the mask image, and the electrode position image.  Specifically, each basis image will" \
  " be read from file \"{BASPREFIX}_b{BASTYPE}_{CHNAME}.fni\", where" \
  " {CHNAME} is the name of the electrode (as specified in the header of" \
  " the input file).  The mask and electrode position images will be read" \
  " from files \"{BASPREFIX}_msk.fni\" and \"{BASPREFIX}_elp.fni\".\n" \
  "\n" \
  "  -skip {NSKIP}\n" \
  "    This optional argument specifies the number of data frames" \
  " to skip from the start of the input EEG file.  If this argument is omitted" \
  " or negative, no frames are skipped.\n" \
  "\n" \
  "  -read {NREAD}\n" \
  "    This optional argument specifies the number of data frames to read" \
  " from the the input EEG file (after skipping the frames specified" \
  " by \"-skip\", if given).  If this argument is omitted," \
  " zero, or negative, the program reads up to the end of the input file.\n" \
  "\n" \
  "  -step {FSTEP}\n" \
  "    This optional argument specifies a subsampling of the input EEG" \
  " data frames.  Namely, an image will be generated only for one out of" \
  " each {FSTEP} consecutive data frames.  If this argument is omitted," \
  " zero, or negative, the program generates one image per data frame.\n" \
  "\n" \
  "  -outPrefix {OUTPREFIX}\n" \
  "    This mandatory argument specifies the common prefix of all output" \
  " images.  Specifically, the image for each data frame will be written" \
  " to the file \"{OUTPREFIX}_b{BASTYPE}_f{KKKKKK}\" where {KKKKKK} is" \
  " the frame index (starting from 0), padded with zeros to six digits.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  neuromat_make_imgbasis(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2013-06-01 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2013-12-02 J. Stolfi: split off {nmeeg_make_imgbasis.c}, substantial rewrite.\n" \
  "  2013-12-02 J. Stolfi: uses {argparser.h}.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " nmeeg_animate_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define stringify(x) strngf(x)
#define strngf(x) #x

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <float_image.h>
#include <float_image_paint.h>
#include <frgb.h>
#include <frgb_path.h>
#include <affirm.h>
#include <jsfile.h>
#include <argparser.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_io.h>
#include <neuromat_eeg_geom.h>
#include <neuromat_image.h>
#include <neuromat_image_png.h>
#include <neuromat_eeg_header.h>
#include <neuromat_eeg_image.h>
#include <neuromat_eeg_image_basis.h>
#include <neuromat_eeg_stats.h>

typedef struct nan_options_t
  { char *basis_prefix; /* Prefix for basis files. */
    char *basis_type;   /* Tag for all basis files. */
    int skip;           /* Data frames to skip. */
    int read;           /* Data frames to read after skipping, or {<=0} for to-end. */
    int step;           /* Animate only one every this many frames. */
    char *outPrefix;    /* Prefix for all output files. */
  } nan_options_t;
  /* Arguments from command line. */
 
float_image_t *nan_read_mask_image(char *prefix);
  /* Reads the mask image from file "{prefix}_msk.fni".  It must be monochromatic,
    with values ranging from 0 (outside nominal scalp) to 1 (inside nominal scalp). */
 
float_image_t *nan_read_electrode_image(char *prefix, int NX, int NY);
  /* Reads  from file "{prefix}_elp.fni" an image that 
    shows the nominal positons of the electrodes, as suitable dots.
    It must be monochromatic, it must have exacty {NX} columns and {NY} rows,
    and its values must range from 0 (transparent) to 1 (opaque) with 
    antialiased edges. */

float_image_t **nan_read_basis_images(char *prefix, char *type, int ne, char *chnames[], int NX, int NY);
  /* Reads a basis {bas[0..ne-1]} of images for interpolation of {ne} electrodes
    whose names are {chname[0..ne-1]}. Each image {bas[ie]} is read from a file called "{prefix}_b{type}_{chname[ie]}".
    It must be monochromatic, it must have exacty {NX} columns and {NY} rows,
    and its values should normally span from 0 (far from electrode {ie})
    to 1 (at the electrode {ie}), but may over- or undershoot that range a little. */ 


float_image_t *nan_make_timeline_bar_image
  ( int NX, 
    int NY, 
    int nt, 
    int nc, 
    double **val, 
    int nm, 
    int ic_mark[], 
    int *bxminP, 
    int *bxmaxP
  );
  /* Creates an RGB image with {NX} columns and {NY} rows containing the
    a timeline bar for {nt} data frames with {nc} channels. Assumes that
    {val[it][ictrig]} is the value of the marker channel in data frame {it}.
    
    For each channel {ic = ic_mark[0..nm-1]}, paints over that timeline bar 
    oneor more rectangles at positions corresponding to the data frames
    where that marker channel is non-zero.
    
    Also returns the range of pixel columns {*bxminP..*bxmaxP} in {0..NX-1} that
    is spanned by the timeline bar. */

void nan_choose_marker_dots_layout
  ( int NXI, 
    int NYI, 
    int NXM, 
    int NYM, 
    int NXB, 
    int NYB, 
    int nm, 
    r2_t ctr_mark[], 
    double *rad_markP
  );
  /* Chooses the the centers
    {ctr_mark[0..nm-1]} and radius {*rad_markP} for the dots that represent
    the {nm} bilevel marker or marker channels in the movie frame images.
    
    Assumes that the main image (showing inerpolated potentials) 
    will have {NXI} columns and {NYI} rows of pixels; below is a bar with size {NXB{ by {NYB};
    and to its left there is a /marker panel/ with width {NXM} and height {NYM}.
    
    If {NXM} and {NYM} are positive, the marker dots are
    arranged vertically in the marker panel, within the rectangle
    {[0 _ NXM] × [0 _ NYI]}. 
    
    If {NXM} or {NYM} is zero, the number {nm} of marker channels must be
    4 or less. In that case the marker dots are placed at the corners of
    the main image. */
        
void nan_paint_marker_channels
  ( float_image_t *frame,
    int nc, 
    double valt[],
    int nm,
    int ic_mark[], 
    double vmax_mark[], 
    r2_t ctr_mark[],
    double rad_mark,
    int style
  );
  /* Paints the marker channels into the given {frame}. Assumes that
    {valt[0..nc-1]} are the sample values of the corresponding data
    frame; that there are {nm} markers to be painted whose channel
    indices are {ic[0..nm-1]}; that the value of marker channel {ic[k]}
    ranges in {[-vmax_mark[k] _ +vmax_mark[k]]}, and it is to be painted
    as a dot with center {ctr_mark[k]} and radius {rad_mark}, with the
    color scale of the given {style}. */

nan_options_t *nan_parse_options(int argc, char **argv);
  /* Parses the command line arguments. */

int main(int argc, char **argv)
  {
    nan_options_t *o = nan_parse_options(argc, argv);
    
    fprintf(stderr, "output files are named %s_f*\n", o->outPrefix);
    fprintf(stderr, "skipping %d frames, reading %d frames\n", o->skip, o->read);
    fprintf(stderr, "animating one every %d frames\n", o->step);
    
    /* Read the EEG header: */
    int nl = 0;
    neuromat_eeg_header_t *h = neuromat_eeg_header_read(stdin, 20, 600, &nl);
    int nc = h->nc; /* Number of channels per data frame (including markers etc.). */
    int ne = h->ne; /* Assumes that the first {ne} channels are electrode potentials. */
    int nm = nc - ne;  /* Number of marker channel dots. */
    fprintf(stderr, "input file has %d channels\n", nc);
    fprintf(stderr, "channels comprise %d electrodes and %d marker channels\n", ne, nm);
    
    /* Read the EEG data frames: */
    /* We must read the frames with step 1, otherwise we may miss trigger pulses. */
    /* But we assume that the electrodes are filtered so that downsampling is aliasing-free.  */
    int nt = 0;
    double **val = neuromat_eeg_data_read(stdin, o->skip, o->read, nc, &nl, &nt);
    fprintf(stderr, "read %d lines, got %d data frames\n", nl, nt);

    /* Compute basic statistics of electrode signals: */
    double vavg[nc], vvar[nc], vmin[nc], vmax[nc];
    double vmin_min, vmax_max, vdev_max;  /* Minimum value among electrode channels. */
    bool_t zeroMean = FALSE;
    neuromat_eeg_stats_per_channel(nt, nc, val, zeroMean, vavg, vvar, vmin, vmax, ne, &vmin_min, &vmax_max, &vdev_max);
    fprintf(stderr, "--- channel statistics ---\n");
    neuromat_eeg_stats_per_channel_print(stderr, nc, h->chname, vavg, vvar, vmin, vmax);
    fprintf(stderr, "range of electrode values = [ %8.5f _ %8.5f ]\n", vmin_min, vmax_max);
    fprintf(stderr, "maximum electrode dev =   %8.5f\n", vdev_max); 
    fprintf(stderr, "\n");

    /* Decide which data frames in file to render, namely {jt_ini..jt_fin} (global) with increment {o->step}: */
    int jt_step = o->step;
    int jt_phase = jt_step/2;  /* First frame in file that is sync'ed to sampling. */
    /* Get {jt_ini} by rounding {o->skip} up to congruent with {jt_phase} mod {jt_step}. */
    int jt_ini = ((o->skip - jt_phase + jt_step - 1)/jt_step)*jt_step + jt_phase;
    assert((jt_ini >= o->skip) && (jt_ini < o->skip + jt_step) && ((jt_ini % jt_step) == jt_phase));
    /* Get {jt_fin} by rounding {o->skip + nt - 1} down to congruent with {jt_phase} mod {jt_step}. */
    int jt_fin = ((o->skip + nt - 1 - jt_phase + jt_step)/jt_step - 1)*jt_step + jt_phase;
    assert((jt_fin < o->skip + nt) && (jt_fin >= o->skip + nt - jt_step) && ((jt_fin % jt_step) == jt_phase));
    demand(jt_fin >= jt_ini, "no output frames in selected data frame range");
    /* Compute number of frames that will be written: */
    int nw = (jt_fin - jt_ini)/jt_step + 1; /* Number of frames to write. */
    fprintf(stderr, "rendering %d frames (%d to %d step %d)\n", nw, jt_ini, jt_fin, jt_step);
      
    /* Read the scalp mask image {msk}, grab main image dimensions {NXI,NYI}: */
    float_image_t *msk = nan_read_mask_image(o->basis_prefix);
    assert(msk->sz[0] == 1);
    int NXI = (int)msk->sz[1];
    int NYI = (int)msk->sz[2];
    
    /* Read the electrode layout mask: */
    float_image_t *elp = nan_read_electrode_image(o->basis_prefix, NXI, NYI);

    /* Read the channel basis images: */
    float_image_t **bas = nan_read_basis_images(o->basis_prefix, o->basis_type, ne, h->chname, NXI, NYI);
    
    /* Allocate the main image (combnation of basis, masked, with electrode overlay):  */
    float_image_t *fim = float_image_new(1, NXI, NYI);

    /* Decide size of marker panel {mpa} at left of main image, if any: */
    int NXM = (nm > 4 ? 20 : 0);   /* Width fo {mpa}: */
    int NYM = (nm > 4 ? NYI : 0);  /* Height of {mpa}: */
    
    /* Decide size of timescale cursor sub-image {bar} below main image, if any: */
    int NXB = (nw > 1 ? NXI : 0); /* Width of {bar} image. */
    int NYB = (nw > 1 ? 20 : 0);  /* Height of {bar}. */
    
    /* Find the max value of marker channels: */
    int ic_mark[nm];
    double vmax_mark[nm];
    int kd;
    for (kd = 0; kd < nm; kd++)
      { int id = ne + kd; /* Index of channel in data frames: */
        ic_mark[kd] = id;
        vmax_mark[kd] = fmax(vmax[id], -vmin[id]);
      }
    
    /* Build image {bar} of timescale bar, with marker channel marks, if any: */
    float_image_t *bar = NULL;
    int bxmin = INT_MIN;   /* Min pixel column of timescale bar in {bar}. */
    int bxmax = INT_MIN;   /* Max pixel column of timescale bar in {bar}. */
    if (nw > 1)
      { assert(NYB > 0);
        bar = nan_make_timeline_bar_image(NXB, NYB, nt, nc, val, nm, ic_mark, &bxmin, &bxmax);
      }
            
    /* Choose position of marker dots in movie frame: */
    r2_t ctr_mark[nm];    /* Centers of marker channel dots. */
    double rad_mark;      /* Radius of marker channel dots. */
    nan_choose_marker_dots_layout(NXI, NYI, NXM, NYM, NXB, NYB, nm, ctr_mark, &rad_mark);

    /* Define the total movie frame dimensions {NX,NY}: */
    int NX = NXM + NXI;
    int NY = NYI + NYB;
    
    /* Paint the frame images:  */
    float_image_t *frame = float_image_new(3, NX, NY); /* The movie frame. */
    double zmax = 1.5; /* Expected max value of interpolated potentials. */
    double coef[ne];
    int nf = 0; /* Frames written. */
    int jt;
    for (jt = jt_ini; jt <= jt_fin; jt += jt_step)
      { int it = jt - o->skip; /* Frame index in {val}. */
        fprintf(stderr, ".");
        
        /* Combine the interp basis images with the channel values as coefs. */
        float_image_fill(fim, 0.0);
        double *valt = val[it];
        int ie;
        for (ie = 0; ie < ne; ie++) 
          { if (nt == 1)
              { /* Single frame: preserve average, scale by max abs electrode value: */ 
                coef[ie] = valt[ie]/fmax(fabs(vmax_max), fabs(vmin_min));
              }
            else
              { /* Multiple frames: mean-shift each channel, divide by max electrode deviation {vdev_max}: */ 
                coef[ie] = (valt[ie] - vavg[ie])/vdev_max;
              }
          }
        neuromat_eeg_image_paint_potentials(ne, coef, bas, msk, 0, fim);

        /* Colorize the image yielding {cim}: */
        int style = 2;
        frgb_t bck = (frgb_t){{ 0.000f, 0.000f, 0.000f }}; /* Background color in mask image. */
        float_image_t *cim = neuromat_image_colorize(fim, msk, zmax, style, &bck);
        
        /* Overlay electrode positions: */
        frgb_t col = (frgb_t){{ 1.000f, 0.900f, 0.000f }}; /* Dot color for channels. */
        neuromat_image_paint_overlay(cim, elp, col.c);
        
        /* Insert main image into movie fram: */
        float_image_assign_rectangle(frame, NXM, NXM+NXI-1, 0, NYI-1, cim, 0, 0);
        float_image_free(cim);
        
        /* Insert time scale and slider dot, if any: */
        if (nt > 1)
          { float_image_assign_rectangle(frame, NXM, NXM+NXB-1, NYI, NYI+NYB-1, bar, 0, 0);
            int it_bar = jt - jt_ini;
            int nt_bar = jt_fin - jt_ini + 1;
            neuromat_image_paint_slider_dot(frame, NXM+bxmin, NXM+bxmax, NYI + 0.5*NYB, (it_bar + 0.5), nt_bar);
          }
        
        /* Paint the marker channels as dots: */
        nan_paint_marker_channels(frame, nc, valt, nm, ic_mark, vmax_mark, ctr_mark, rad_mark, style);

        char *tag = NULL;
        asprintf(&tag, "f%06d", jt);
        neuromat_image_png_write(o->outPrefix, tag, frame, 0.0, 1.0);
        free(tag);
        nf++;
        if ((nf % 50) == 0) { fprintf(stderr, "\n"); }
      }
    if ((nf % 50) != 0) { fprintf(stderr, "\n"); }
    assert(nf == nw);
 

    float_image_free(frame);
    if (bar != NULL) { float_image_free(bar); }
    float_image_free(fim);
    float_image_free(msk);
    float_image_free(elp);
    { int ie; for (ie = 0; ie < ne; ie++) { float_image_free(bas[ie]); }  }
    free(bas);
    { int it; for (it = 0; it < nt; it++) { free(val[it]); } }
    free(val);
    return 0;
  }

float_image_t *nan_read_mask_image(char *prefix)
  { char *fname = NULL;
    asprintf(&fname, "%s_msk.fni", prefix);
    FILE *rd = open_read(fname, TRUE);
    float_image_t *msk = float_image_read(rd);
    fclose(rd);
    demand(msk->sz[0] == 1, "mask image must have a single channel");
    free(fname);
    return msk;
  }

float_image_t *nan_read_electrode_image(char *prefix, int NX, int NY)
  { char *fname = NULL;
    asprintf(&fname, "%s_elp.fni", prefix);
    FILE *rd = open_read(fname, TRUE);
    float_image_t *elp = float_image_read(rd);
    fclose(rd);
    demand(elp->sz[0] == 1, "electrode position image must have a single channel");
    demand(elp->sz[1] == NX, "electrode position image has wrong width");
    demand(elp->sz[2] == NY, "electrode position image has wrong height");
    free(fname);
    return elp;
  }

float_image_t **nan_read_basis_images(char *prefix, char *type, int ne, char *chname[], int NX, int NY)
  {
    float_image_t **bas = notnull(malloc(ne*sizeof(float_image_t*)), "no mem");
    int ie;
    for (ie = 0; ie < ne; ie++)
      { char *fname = NULL;
        asprintf(&fname, "%s_b%s_%s.fni", prefix, type, chname[ie]);
        FILE *rd = open_read(fname, TRUE);
        float_image_t *bsi = float_image_read(rd);
        fclose(rd);
        demand(bsi->sz[0] == 1, "electrode position image must have a single channel");
        demand(bsi->sz[1] == NX, "electrode position image has wrong width");
        demand(bsi->sz[2] == NY, "electrode position image has wrong height");
        bas[ie] = bsi;
        free(fname);
      }
    return bas;
  }

void nan_paint_marker_channels
  ( float_image_t *frame,
    int nc, 
    double valt[],
    int nm,
    int ic_mark[], 
    double vmax_mark[], 
    r2_t ctr_mark[],
    double rad_mark,
    int style
  )
  { int km;
    for (km = 0; km < nm; km++)
      { int ic = ic_mark[km]; /* Index of channel in data frames: */
        demand((ic >= 0) && (ic < nc), "invalid marker channel");
        /* Convert marker value to {[-1 _ +1]}: */
        float zk = (float)(valt[ic]/vmax_mark[km]); /* Marker value scaled to {[-1..+1]} */
        frgb_t rgbk = frgb_path_map_signed(zk, 1, style);
        int c;
        for (c = 0; c < 3; c++)
          { r2_t *ctr = &(ctr_mark[km]);
            float_image_paint_dot(frame, c, ctr->c[0], ctr->c[1], rad_mark, 0, TRUE, FALSE, rgbk.c[c], NAN, 4);
          }
      }
  }

float_image_t *nan_make_timeline_bar_image
  ( int NX, 
    int NY, 
    int nt, 
    int nc, 
    double **val, 
    int nm, 
    int ic_mark[], 
    int *bxminP, 
    int *bxmaxP
  )
  {
    demand(NY % 4 == 0, "marker image height must be multiple of 4");
    demand(NY >= 8, "marker image height too small");
    demand(NX > 3*NY/2, "marker image aspect too small");
    
    /* Build an image with the timeline bar, and get its column span {*bxminP..*bxmaxP}: */
    float_image_t *img = neuromat_image_make_slider_bar(3, NX, NY, bxminP, bxmaxP);
    
    /* Choose row span {pymin..pymax} of the fattest markers: */
    int pysize = NY - 4; /* Leave 2 pixels above and below. */
    int pymin = (NY - pysize)/2;
    int pymax = NY - 1 - pymin;
    
    /* Choose row span {qymin..qymax} of the skinniest markers: */
    int qysize = 4; /* Min 4 pixels. */
    assert(qysize <= pysize);
    int qymin = (NY - qysize)/2;
    int qymax = NY - 1 - qymin;
    
    /* Int how many sizes do we have? */
    int nsizes = pymax - qymax + 1;
    
    int km;
    for (km = 0; km < nm; km++)
      { /* Paint the marker marks: */
        int ic = ic_mark[km];
        demand((ic >= 0) && (ic < nc), "invalid marker channel");
        int mymin = qymin + (km % nsizes);
        int mymax = qymax - (km % nsizes);
        neuromat_eeg_image_paint_marker_tics(img, *bxminP, *bxmaxP, mymin, mymax, nt, nc, val, ic);
      }
    return img;
  }
    
void nan_choose_marker_dots_layout
  ( int NXI, 
    int NYI,
    int NXM, 
    int NYM, 
    int NXB, 
    int NYB, 
    int nm, 
    r2_t ctr_mark[], 
    double *rad_markP
  )
  {
    /* Column of marker lights is unneeded if less than 4 lights: */
    bool_t marker_bar = ((NXM > 0) && (NYM > 0)); /* True iff uses a marker dot column on left side. */ 
    if (! marker_bar)
      { NXM = 0; NYM = 0;
        demand(nm <= 4, "too many markers for corner placement");
      }
    
    if (marker_bar)
      { /* Put a column of dots in a separate bar at left. */
        /* Get center X {ctx} and Y step {dty}: */
        double ctx = 0.5*NXM; /* Center X of dots. */
        double dty = ((double)NYM)/nm; /* Y spacing of dot centers. */
        /* Position the dots: */
        int km;
        for (km = 0; km < nm; km++) { ctr_mark[km] = (r2_t){{ ctx, NYM - (km + 0.5)*dty }}; }
        (*rad_markP) = 0.45*fmin((double)NXM, dty);
      }
    else
      { /* Put marker dots at corners of main image, in reading order. */
        /* Assumes that scalp diagram fits in the ellipse inscreibed in main image domain. */
        /* Find the largest proportional rectangle that will fit between ellipse and corner: */
        double hx = 0.5*NXI;
        double hy = 0.5*NYI;
        double hmin = fmin(hx, hy);
        
        /* Now assume that the main image is a square {S} with side {2*hmin}. */
        double wmin = (M_SQRT2 - 1.0)*hmin; /* Dist from corner of {S} to inscribed circle. */
        double rad = wmin/(1 + M_SQRT2);    /* Radius of circle that fits snugly in that corner. */
        /* Place the dot centers: */
        int idx, idy;
        for (idy = 0; idy < 2; idy++)
          { double cty = (idy == 0 ? NYB + NYI - rad : NYB + rad);
            for (idx = 0; idx < 2; idx++)
              { double ctx = (idx == 0 ? NXM + rad : NXM + NXI - rad);
                int km = 2*idy + idx;
                if (km < nm) { ctr_mark[km] = (r2_t){{ ctx, cty }}; }
              }
          }
        (*rad_markP) = 0.90*rad;
      }
  }

nan_options_t *nan_parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    /* Allocate the command line argument record: */
    nan_options_t *o = notnull(malloc(sizeof(nan_options_t)), "no mem");

    /* Parse keyword parameters: */

    argparser_get_keyword(pp, "-basis");
    o->basis_prefix = argparser_get_next_non_keyword(pp);
    o->basis_type = argparser_get_next_non_keyword(pp);
    

    if (argparser_keyword_present(pp, "-skip"))
      { o->skip = (int)argparser_get_next_int(pp, 0, INT_MAX); }
    else
      { o->skip = 0; }
      
    if (argparser_keyword_present(pp, "-read"))
      { o->read = (int)argparser_get_next_int(pp, INT_MIN, INT_MAX);
        if (o->read <= 0) { o->read = INT_MIN; }
      }
    else
      { o->read = INT_MIN; }
      
    if (argparser_keyword_present(pp, "-step"))
      { o->step = (int)argparser_get_next_int(pp, 1, INT_MAX); }
    else
      { o->step = 1; }
      
    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next_non_keyword(pp);

    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }
