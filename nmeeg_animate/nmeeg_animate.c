#define PROG_NAME "nmeeg_animate"
#define PROG_DESC "Produces a colored animation of an EEG dataset."
#define PROG_VERS "2021-08-24"

#define nmeeg_animate_C_COPYRIGHT \
  "Copyright © 2013 by the State University of Campinas (UNICAMP)"
/* Last edited on 2023-10-22 11:29:32 by stolfi */

/* !!! Colorize marker channels discontinuously around zero. !!! */

/* !!! There should be separate {skip,read,step} for each dataset. !!! */ 

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -capType {CAP_TYPE} \\\n" \
  "    -basisDir {BAS_DIR} \\\n" \
  "    -maskDir {MSK_DIR} \\\n" \
  "    [ -skip {NSKIP} ] \\\n" \
  "    [ -read {NREAD} ] \\\n" \
  "    [ -step {FSTEP} ] \\\n" \
  "    -inDir {IN_PREFIX} \\\n" \
  "    { -setName {SET_NAME} }.. \\\n" \
  "    { -marker {MK_CH} {MK_VAL} {MK_R} {MK_G} {MK_B} }.. \\\n" \
  "    -outDir {OUT_DIR} \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program reads one or more EEG datasets" \
  " and writes an inage for each frame of each dataset, showing the electrode intensities" \
  " over the idealized 2D scalp as interpolated from the electrode potential, overlaying an idealized diagram of" \
  " the scalp.  These images can be converted into an animation of the dataset.\n" \
  "\n" \
  "  The program assumes that every data line (/data frame/) in each input dataset contains" \
  " simultaneous samples from some fixed number {NC} of data channels, of which" \
  " the first {NE} are electrode potentials, while the rest are phase indicators (/markers/)" \
  " or other unspecified signals.\n" \
  "\n" \
  "  The program first reads a /mask image/ {MSK}, a monochromatic image that defines the part of the" \
  " image that corresponds to the region of the scalp that is to be" \
  " displayed; an /electrode mask/ {ELP} that shows the position of each electrode over that region; and" \
  " an /image basis/, a set of {NE} single-channel image {BIM[0..NE-1]} that define" \
  " how each electrode value is to be rendered.\n" \
  "\n" \
  "  For each input frame {VI} in each input dataset, the program combines the images {BIM[0..NE-1]} pixel by pixel" \
  " into a single monochromatic image {IMG}, using the electrode values as coefficients.  That is,\n" \
  "\n" \
  "   {IMG = VI[0]*BIM[0] + VI[1]*BIM[1] + ··· + VI[NE-1]*BIM[NE-1]}\n" \
  "\n" \
  " The program will then convert each pixel" \
  " to an RGB color value according to a built-in false-color scale, overlay the electrode" \
  " position image {ELP} to it, insert the {MSK} image as the alpha channel, and output {IMG} to a" \
  " separate file.\n" \
  "\n" \
  "  If marker channels are specified (see the \"-marker\" option), each frame image {IMG} will" \
  " also have a set of /time tracks/ appended at the bottom, showing the" \
  " position of the frame in the dataset, and the ranges in that dataset where specified marker" \
  " channels are on. The states of those marker channels is also shown by color-coded dots" \
  " to its left of the time tracks.  Note that, when markers are specified, the final" \
  " frane image will be taller than the basis and mask images.\n" \
  "\n" \
  "FILE FORMATS\n" \
  "  The basis, mask, and electrode position files must be float-valued images" \
  " as written by {float_image_write}. All images must have a single channel and the same" \
  " size (rows and columns) in pixels.\n" \
  "\n" \
  "  The basis images may have positive and negative values.\n" \
  "\n" \
  "  The mask image samples must be 0 in the background, 1 inside the" \
  " idealized head.  The electrode position image must be 0 for" \
  " transparent, 1 for opaque.  Both should be antialiased (between" \
  " 0 and 1) along the boundary.\n" \
  "\n" \
  "  Output image files will be in PNG format.  If marker channels are" \
  " specified (see the \"-marker\" option), the frame height will be more than  \n" \
  "\n" \
  "OPTIONS\n" \
  "  -capType {CAP_TYPE}\n" \
  "    This mandatory argument specifies the cap type, which determines" \
  " the count, names, and positions of the electrodes.  The possible" \
  " values of {CAP_TYPE} are defined" \
  " by {neuromat_eeg_get_channel_names} in {neuromat_eeg.h}, and" \
  " include \"R20\", \"R128\", \"R129\", and \"FN3\".\n" \
  "\n" \
  "  -basisDir {BAS_DIR}\n" \
  "    This mandatory argument specifies the directory where the basis" \
  " images are to be found.  Specifically, each basis image will" \
  " be read from file \"{BAS_DIR}/{CHNAME}.fni\", where" \
  " {CHNAME} is the name of the electrode (as specified in the header of" \
  " the input file).\n" \
  "\n" \
  "  -maskDir {MSK_DIR}\n" \
  "    This mandatory argument specifies the directory where the scalp mask image and" \
  " the electrode positions overlay image are to be found.  Specifically, they will be read" \
  " from files \"{MSK_DIR}/msk.fni\" and \"{MSK_DIR}/elp.fni\".\n" \
  "\n" \
  "  -skip {NSKIP}\n" \
  "    This optional argument specifies the number of data frames" \
  " to skip from the start of each input EEG file.  If this argument is omitted" \
  " or negative, no frames are skipped.\n" \
  "\n" \
  "  -read {NREAD}\n" \
  "    This optional argument specifies the number of data frames to read" \
  " from each input EEG file (after skipping the frames specified" \
  " by \"-skip\", if given).  If this argument is omitted," \
  " zero, or negative, the program reads up to the end of the input file.\n" \
  "\n" \
  "  -step {FSTEP}\n" \
  "    This optional argument specifies a subsampling of the input EEG" \
  " data frames.  Namely, an image will be generated only for one out of" \
  " each {FSTEP} consecutive data frames.  If this argument is omitted," \
  " zero, or negative, the program generates one image per data frame.\n" \
  "\n" \
  "  -inDir {IN_DIR}\n" \
  "    This mandatory argument specifies the directory where the input EEG dataset files are found." \
  " \n" \
  "  -outDir {OUT_DIR}\n" \
  "    This mandatory argument specifies the directory wher the output" \
  " frame images will be written to.\n" \
  "  -setName {SET_NAME}\n" \
  "    This mandatory argument specifies the name of the input dataset (without extension)." \
  " The EEG data will be read from \"{IN_DIR}/{SET_NAME}.txt\" and the frame images will be" \
  " written out to \"{OUT_DIR}/{SET_NAME}/f{KKKKKK}.fni\", where {KKKKKK} is" \
  " the frame index (starting from 0), padded with zeros to six digits.  This argument may be" \
  " repeated to process more than one dataset with the same basis images.\n" \
  "  -marker {MK_CH} {MK_VAL} {MK_R} {MK_G} {MK_B}\n" \
  "    This optional argument specifies a marker channel to be indicated in the output" \
  " frames. It may be repeated zero or more times; for each time, a gray /time track/ line" \
  " will be drawn onto the frame images, below the interpolated potentials" \
  " image." \
  " \n" \
  "    If {MK_CH} is in the range {0..NC-1}, the marker value {V} will be" \
  " taken from the data channel with that index.  The frame ranges where those" \
  " values are nonzero will be highlighted over the time track with RGB" \
  " color {({MK_R},{MK_G},{MK_B})}." \
  " \n" \
  "    In addition, a dot will be painted with" \
  " that color times {V/MK_VAL} if that ratio is positive, and the" \
  " complementary color times {-V/MK_VAL} if that ratio is negative.  In particular, the" \
  " dot will be black if {V} is zero.  These dots will be painted in to the" \
  " left of the corresponding time tracks." \
  " \n" \
  "    If {MK_CH} is outside the" \
  " range {0..NC-1}, its value {V} will be assumed to be zero, so the time track will have no highlight" \
  " and its dot will always be black.\n" \
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
  "  2017-09-27 J. Stolfi: substantial rewrite of the timeline bar.\n" \
  "  2021-08-23 J. Stolfi: changed \"-basis\" option to string-valued {BAS_NAME}.\n" \
  "  2021-08-28 J. Stolfi: changed \"-basPrefix\" \"-inPrefix\" \"-outPrefix\" to \"-xxxDir\".\n" \
  "  2021-08-28 J. Stolfi: separated \"-maskPrefix\" from \"-basPrefix\".\n" \
  "  2021-08-28 J. Stolfi: Changed output file names to one directory per dataset.\n" \
  "  2021-08-28 J. Stolfi: Added the \"-marker\" option.\n" \
  "  2021-08-31 J. Stolfi: changed \"-electrodes\" from count to list of names.\n" \
  "  2023-10-22 J. Stolfi: replaced \"-electrodes\" by \"-capType\".\n" \
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
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <float_image.h>
#include <float_image_paint.h>
#include <float_image_overlay.h>
#include <frgb.h>
#include <frgb_path.h>
#include <affirm.h>
#include <vec.h>
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
#include <neuromat_eeg_channel_stats.h>

typedef neuromat_eeg_marker_spec_t nan_mkop_t;
  /* Data parsed from a "-marker" option. */
  
vec_typedef(nan_mkop_vec_t,nan_mkop_vec,nan_mkop_t);

typedef struct nan_options_t
  { char *capType;              /* Type of EEG cap, like "R20", "R128" etc. */
    char *basisDir;          /* Dir for basis image files. */
    char *maskDir;           /* Dir for mask and electrode overlay images. */
    int32_t skip;            /* Data frames to skip. */
    int32_t read;            /* Data frames to read after skipping, or {<=0} for to-end. */
    int32_t step;            /* Animate only one every this many frames. */
    char *inDir;             /* Dir for all input EEG dataset files. */
    string_vec_t setName;    /* Names of EEG datasets to animate. */
    nan_mkop_vec_t marker;   /* Marker channel specs. */
    char *outDir;            /* Dir for all output files. */
  } nan_options_t; 
  /* Arguments from command line. */
 
/* NOTE: Float images conventionally have the Y axis pointing UP with row 0 at BOTTOM. */

float_image_t *nan_read_mask_image(char *maskDir);
  /* Reads the mask image from file "{maskDir}/msk.fni".  It must be monochromatic,
    with values ranging from 0 (outside the relevant portion of the scalp) to
    1 (inside). */
 
float_image_t *nan_read_electrode_image(char *maskDir, int32_t NX, int32_t NY);
  /* Reads  from file "{maskDir}/elp.fni" an image that 
    shows the nominal positons of the electrodes, as suitable dots.
    It must be monochromatic, it must have exacty {NX} columns and {NY} rows,
    and its values must range from 0 (transparent) to 1 (opaque) with 
    antialiased edges. */

float_image_t **nan_read_basis_images(char *basisDir, int32_t ne, char *chname[], int32_t NX, int32_t NY);
  /* Reads a basis {bas[0..ne-1]} of images for interpolation of {ne}
    electrodes whose names are {chname[0..ne-1]}. Each image {bas[ie]}
    is read from a file called "{basisDir}/{chname[ie]}". It must
    be monochromatic, it must have exacty {NX} columns and {NY} rows,
    and its values should normally span from 0 (far from electrode {ie})
    to 1 (at the electrode {ie}), but may over- or undershoot that range
    a little. */ 

void nan_animate_dataset
  ( char *inDir, 
    char *setName, 
    int32_t skip, 
    int32_t read, 
    int32_t step, 
    int32_t ne,
    float_image_t *msk, 
    float_image_t *elp, 
    float_image_t **bas, 
    nan_mkop_vec_t *marker,
    char *outDir
  );
  /* Reads a dataset from file "{inDir}/{setName}.txt" and writes each selected frame
    as a PNG image to "{outDir}/{setName}/f{KKKKKK}.png". 
    
    The program skips the first {skip} frames, then reads {read} frames, and processes only
    the first of every {step} frames read. The dataset had better have {ne} electrodes.

    The output images are created by combning the images {msk}, {elp} and {bas[0..ne-1]} and 
    the marker time tracks and dots specified by {marker},  as described in the {PROG_INFO} string  */

void nan_choose_channel_scaling_parms
  ( int32_t nt,         /* Number of data frames. */
    int32_t nc,         /* Number of data channels (including markers). */
    double **val,       /* Sample values per frame and channel. */
    char **chname,     /* Channel names. */
    int32_t ne,         /* Number of electrodes. */
    double vscale[]     /* (OUT) value to be divided into each channel value. */
  );
  /* Chooses suitable scaling parameters for the channel values.
    Assumes the channel samples are {val[it][ic]} for {it} in
    {0..nt-1} and {ic} in [0..nc-1}. Assumes that the first {ne}
    channels are electrode potentials, and the rest are marker
    channels.

    Returns in {vscale[ic]}, for {ic} in {0..nc-1}, the amount to
    divide into each sample of channel {ic} before colorizing it.
    Also prints channel statistics to {stderr}.

    The scaling factors {vscale[ic]} are chosen differently
    depending on whether channel {ic} is an electrode potential (in
    {0..ne-1}) or a marker (in {ne..nc-1}). They depend on the RMS
    value {vrms[ic]} of the channel samples {val[0..nt-1][ic]}.

    For an electrode channel {ie} in {0..ne-1}: if there is only one
    frame, then {vscale[ie]} is the maximum absolute value of all
    electrode potential samples, otherwise it is the maximum of
    {vrms[0..ne-1]}.

    For a marker channel {ic} in {ne..nc-1}: {vscale[ic]} is the max
    absolute sample value of that channel. */
        
void nan_make_main_frame_image
  ( int32_t nc,           /* Number of channels per frame. */
    double valt[],        /* Channel sample values (µV). */
    int32_t ne,           /* Number of electrodes. */
    double vscale[],      /* Scaling denominator for each channel. */
    float_image_t *bas[], /* Basis images (signed, monochromatic). */
    float_image_t *msk,   /* Mask image with nominal outline of scalp. */
    float_image_t *elp,   /* Overlay image showing electrode positions. */
    float_image_t *fimg_mono,   /* WORK: Temporary monochrome frame image. */
    float_image_t *fimg_rgba    /* OUT: Colorized frame image. */
  );
  /* Fills the image {fimg_rgba} with a colorized visualization of
    the potentials {valt[0..ne-1]}, mapped to an idealized scalp
    shape by the image basis {bas[0..ne-1]}.

    The procedure assumes that {val[0..nc-1]} are the channel
    samples of one frame, of which the first {ne} are electrode potentials.
    Each sample {val[ie]}, for {ie} in {0..ne-1},
    is divided by {vscale[ie]}.  The shifted and scaled samples
    are then used as coefficients to combine the basis
    images {bas[0..ne-1]}.  The result is then colorized,
    clipped by the mask image {msk}, and the image {elp}
    is overlaid on it. 

    The result is returned in the image {fimg_rgba} that should have 3
    channels.  The single-channel image {fimg_mono} is used as a 
    scratch area.  All images must have the same number of 
    columns and rows. */

nan_options_t *nan_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments. */

int32_t main(int32_t argc, char **argv)
  {
    nan_options_t *o = nan_parse_options(argc, argv);
      
    int32_t ne = -1;
    char **chname = NULL;
    neuromat_eeg_get_channel_names(o->capType, 0, NULL, &ne, &chname);

    /* Read the scalp mask image {msk}, grab main image dimensions {NXI,NYI}: */
    float_image_t *msk = nan_read_mask_image(o->maskDir);
    assert(msk->sz[0] == 1);
    int32_t NXI = (int32_t)msk->sz[1];
    int32_t NYI = (int32_t)msk->sz[2];
    
    /* Read the electrode layout mask: */
    float_image_t *elp = nan_read_electrode_image(o->maskDir, NXI, NYI);
 
    /* Read the channel basis images: */
    float_image_t **bas = nan_read_basis_images(o->basisDir, ne, chname, NXI, NYI);
    
    /* Process the dataset: */
    int32_t nd = o->setName.ne; /* Number of datasets to process. */
    for (int32_t id = 0; id < nd; id++)
      { char *setName = o->setName.e[id];
        nan_animate_dataset(o->inDir, setName, o->skip, o->read, o->step, ne, msk, elp, bas, &(o->marker), o->outDir);
      }

    float_image_free(msk);
    float_image_free(elp);
    { int32_t ie; for (ie = 0; ie < ne; ie++) { float_image_free(bas[ie]); }  }
    free(bas);
    
    return 0;
  }
 
void nan_animate_dataset
  ( char *inDir, 
    char *setName, 
    int32_t skip, 
    int32_t read, 
    int32_t step, 
    int32_t ne,
    float_image_t *msk, 
    float_image_t *elp, 
    float_image_t **bas, 
    nan_mkop_vec_t *marker, 
    char *outDir
  )
  {
    assert(msk->sz[0] == 1);
    int32_t NXI = (int32_t)msk->sz[1];
    int32_t NYI = (int32_t)msk->sz[2];
    
    fprintf(stderr, "skipping %d frames, reading %d frames\n", skip, read);
    fprintf(stderr, "animating one every %d frames\n", step);

    /* Open the file: */
    char *inFile = NULL;
    asprintf(&inFile, "%s/%s.txt", inDir, setName);
    fprintf(stderr, "processing the dataset %s\n", inFile);
    FILE *rd = open_read(inFile, TRUE);
    free(inFile);
    
    /* Read the EEG header: */
    int32_t nl = 0;
    neuromat_eeg_header_t *h = neuromat_eeg_header_read(rd, 20, 600, &nl);
    int32_t nc = h->nc; /* Number of channels per data frame (including markers etc.). */
    demand(h->ne == ne, "wrong number of electrodes in input dataset");
    fprintf(stderr, "input file has %d channels\n", nc); 
    fprintf(stderr, "the first %d channels are assumed to be electrode potentials\n", ne);
    
    /* Read the EEG data frames: */
    /* We must read the frames with step 1, otherwise we may miss trigger pulses. */
    /* But we assume that the electrodes are filtered so that downsampling is aliasing-free.  */
    int32_t nt = 0;
    double **val = neuromat_eeg_data_read(rd, skip, read, nc, &nl, &nt);
    fprintf(stderr, "read %d lines, got %d data frames\n", nl, nt);
    fclose(rd);

    /* Choose value scaling params for electrode potentials: */
    double vscale[ne]; /* Scale to use for electrode potentials. */
    nan_choose_channel_scaling_parms(nt, nc, val, h->chname, ne, vscale);
    
    /* Decide which data frames in file to render, {jt_ini..jt_fin} (global) with increment {step}: */
    int32_t jt_step = step;
    int32_t jt_phase = jt_step/2;  /* First frame in file that is sync'ed to sampling. */
    /* Get {jt_ini} by rounding {skip} up to congruent with {jt_phase} mod {jt_step}. */
    int32_t jt_ini = ((skip - jt_phase + jt_step - 1)/jt_step)*jt_step + jt_phase;
    assert((jt_ini >= skip) && (jt_ini < skip + jt_step) && ((jt_ini % jt_step) == jt_phase));
    /* Get {jt_fin} by rounding {skip + nt - 1} down to congruent with {jt_phase} mod {jt_step}. */
    int32_t jt_fin = ((skip + nt - 1 - jt_phase + jt_step)/jt_step - 1)*jt_step + jt_phase;
    assert((jt_fin < skip + nt) && (jt_fin >= skip + nt - jt_step) && ((jt_fin % jt_step) == jt_phase));
    demand(jt_fin >= jt_ini, "no output frames in selected data frame range");
    /* Compute number of frames that will be written: */
    int32_t nw = (jt_fin - jt_ini)/jt_step + 1; /* Number of frames to write. */
    fprintf(stderr, "rendering %d frames (%d to %d step %d)\n", nw, jt_ini, jt_fin, jt_step);
    
    /* Allocate the main image (combnation of basis, masked, with electrode overlay):  */
    float_image_t *fimg_mono = float_image_new(1, NXI, NYI); /* Working single-channel numeric image. */
    float_image_t *fimg_rgba = float_image_new(4, NXI, NYI); /* Colorized version of {fimg_mono}. */

    /* Create the time track and dots image for the marker channels: */
    int32_t nm = marker->ne;  /* Number of marker channel dots. */
    fprintf(stderr, "assuming %d marker channels\n", nm);
    
    /* Define the line width for time tracks, slider: */
    int32_t hwB = (int32_t)imax(1, (int32_t)floor(NXI/250.0 + 0.5)); /* Line half-width in pixels. */
    
    float_image_t *bar = NULL; /* Image with time tracks and dots, or NULL if none. */
    r2_t mkdots_ctr[nm];        /* Centers of marker channel dots. */
    double mkdots_rad = -1;     /* Radius of marker channel dots. */
    int32_t track_xlo = -1;     /* Left offset of timeline tracks {bar}. */
    int32_t track_xsz = -1;     /* Dimensions of timeline in {bar}, if any. */
    int32_t track_y[nm];        /* Y offset to centerline of each track. */
    int32_t slider_hh;          /* Y extent slider marker, from timeline. */
    if (nm > 0)
      { /* Decide size of marker time track and dots sub-image {bar} below main image, if any: */
        int32_t mkdots_xctr = -1;      /* X position of dot centers. */
        bar = neuromat_eeg_image_make_time_tracks
          ( nt, nc, val, nm, marker->e, hwB, NXI, &mkdots_xctr, &mkdots_rad, &track_xlo, &track_xsz, track_y, &slider_hh );
        assert(bar->sz[0] == 4);
        assert(bar->sz[1] == NXI);
        /* Define dot centers: */
        for (int32_t im = 0; im < nm; im++) { mkdots_ctr[im] = (r2_t){{ mkdots_xctr, track_y[im] }}; }
      }

    /* Define the total movie frame dimensions {NX,NY}: */
    int32_t NYB = (bar == NULL ? 0 : (int32_t)bar->sz[2]);
    int32_t NX = NXI;
    int32_t NY = NYI + NYB;
    
    char *frameDir = NULL;
    asprintf(&frameDir, "%s/%s", outDir, setName);
    fprintf(stderr, "output files are named %s/f{KKKKKK}.png\n", frameDir);
    
    /* Paint the frame images:  */
    float_image_t *frame = float_image_new(4, NX, NY); /* The movie frame. */
    double tlo = 0.0; /* Time value for start of timeline. */
    double thi = (double)nt; /* Time value for end of timeline. */
    frgb_t fslider = (frgb_t){{ 1.000f, 0.800f, 0.000f }}; /* Color of slider. */
    int32_t nf = 0; /* Frames written. */
    for (int32_t jt = jt_ini; jt <= jt_fin; jt += jt_step)
      { int32_t it = jt - skip; /* Frame index in {val}. */
        fprintf(stderr, ".");
        double *valt = val[it];
        nan_make_main_frame_image(nc, valt, ne, vscale, bas, msk, elp, fimg_mono, fimg_rgba);
        
        /* Insert main image into movie frame: */
        float_image_assign_rectangle(frame, 0, NXI-1, NYB, NYB+NYI-1, fimg_rgba, 0, 0);
        
        /* Insert time track panel, if any: */
        if (bar != NULL)
          { /* Copy timeline bar image into full image: */
            float_image_assign_rectangle(frame, 0, NXI-1, 0, NYB-1, bar, 0, 0);
            /* Paint the slider over it: */
            assert(nm > 0);
            double t = it + 0.5; /* Nominal time of frame, in framesteps. */
            neuromat_image_paint_slider
              ( frame, hwB, slider_hh - hwB, tlo, t, thi, track_xlo, track_xsz, track_y[0] - hwB, track_y[nm-1] + hwB, fslider );
            /* Paint the marker channels as dots: */
            neuromat_eeg_image_paint_marker_dots
              ( frame, nc, valt, nm, marker->e, mkdots_ctr, mkdots_rad );
          }

        char *frameName = NULL;
        asprintf(&frameName, "f%06d", jt);
        neuromat_image_png_write( frameDir, frameName, frame, 0.0, 1.0, NAN );
        free(frameName);
        nf++;
        if ((nf % 50) == 0) { fprintf(stderr, "\n"); }
      }
    if ((nf % 50) != 0) { fprintf(stderr, "\n"); }
    assert(nf == nw);
 
    free(frameDir);
    float_image_free(frame);
    if (bar != NULL) { float_image_free(bar); }
    float_image_free(fimg_mono);
    float_image_free(fimg_rgba);
    { int32_t it; for (it = 0; it < nt; it++) { free(val[it]); } }
    free(val);
  }

float_image_t *nan_read_mask_image(char *maskDir)
  { char *fname = NULL;
    asprintf(&fname, "%s/msk.fni", maskDir);
    FILE *rd = open_read(fname, TRUE);
    float_image_t *msk = float_image_read(rd);
    fclose(rd);
    demand(msk->sz[0] == 1, "mask image must have a single channel");
    free(fname);
    return msk;
  }

float_image_t *nan_read_electrode_image(char *maskDir, int32_t NX, int32_t NY)
  { char *fname = NULL;
    asprintf(&fname, "%s/elp.fni", maskDir);
    FILE *rd = open_read(fname, TRUE);
    float_image_t *elp = float_image_read(rd);
    fclose(rd);
    demand(elp->sz[0] == 4, "electrode position image should be RGBA");
    demand(elp->sz[1] == NX, "electrode position image has wrong width");
    demand(elp->sz[2] == NY, "electrode position image has wrong height");
    free(fname);
    return elp;
  }

float_image_t **nan_read_basis_images(char *basisDir, int32_t ne, char *chname[], int32_t NX, int32_t NY)
  {
    float_image_t **bas = notnull(malloc(ne*sizeof(float_image_t*)), "no mem");
    int32_t ie;
    for (ie = 0; ie < ne; ie++)
      { char *fname = NULL;
        asprintf(&fname, "%s/%s.fni", basisDir, chname[ie]);
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

void nan_choose_channel_scaling_parms
  ( int32_t nt,             /* Number of data frames. */
    int32_t nc,             /* Number of data channels (including markers). */
    double **val,       /* Sample values per frame and channel. */
    char **chname,     /* Channel names. */
    int32_t ne,             /* Number of electrodes. */
    double vscale[]     /* (OUT) value to be divided into each channel value. */
  )
  {
    /* Compute basic statistics of electrode signals: */
    neuromat_eeg_channel_stats_t *st = neuromat_eeg_channel_stats_new(nc);
    neuromat_eeg_channel_stats_t *stg = neuromat_eeg_channel_stats_new(1);
    double eps = 0.01; /* Assumed uncertainty of measurement (µV). */
    neuromat_eeg_channel_stats_gather_all(nt, nc, val, NULL, eps, st, ne, stg);
    fprintf(stderr, "  --- channel statistics ---\n");
    neuromat_eeg_channel_stats_print_all(stderr, 2, nc, chname, FALSE, st, ne, stg);
    fprintf(stderr, "\n");

    /* Select scaling factors for electrode potentials: */
    double escale;
    if (nt == 1)
      { /* Single frame: scale by max abs electrode value: */ 
        escale = fmax(stg->max, -stg->min);
      }
    else
      { /* Multiple frames: Scale by max rms electrode value: */ 
        escale = stg->rms;
      }
    for (int32_t ie = 0; ie < ne; ie++) 
      { vscale[ie] = escale; }

    /* Select scaling factors for marker channels: */
    for (int32_t ic = ne; ic < nc; ic++)
      { vscale[ic] = fmax(st[ic].max, -st[ic].min); }
      
    /* Sanity check: */
    for (int32_t ic = 0; ic < nc; ic++) 
      { if ((isnan(vscale[ic])) || (vscale[ic] <= 0)) 
          { fprintf(stderr, "** wscale[%d] = %24.16f\n", ic, vscale[ic]); demand(FALSE, "aborted"); }
      }

    free(st);
    free(stg);
  }

void nan_make_main_frame_image
  ( int32_t nc,           /* Number of channels per frame. */
    double valt[],        /* Channel sample values (µV). */
    int32_t ne,           /* Number of electrodes. */
    double vscale[],      /* Scaling denominator for each channel. */
    float_image_t *bas[], /* Basis images (signed, monochromatic). */
    float_image_t *msk,   /* Mask image with nominal outline of scalp. */
    float_image_t *elp,   /* Overlay image showing electrode positions. */
    float_image_t *fimg_mono,   /* WORK: Temporary monochrome frame image. */
    float_image_t *fimg_rgba    /* OUT: Colorized frame image. */
  )
  {
    double zmax = 1.5; /* Expected max value of interpolated potentials. */

    /* Combine the interp basis images with the channel values as coefs. */
    float_image_fill(fimg_mono, 0.0);
    double coef[ne];
    for (int32_t ie = 0; ie < ne; ie++) { coef[ie] = valt[ie]/vscale[ie]; }
    neuromat_eeg_image_compute_pot_field(ne, coef, bas, msk, fimg_mono);

    /* Colorize the image yielding {fimg_rgba}: */
    int32_t style = 2;
    neuromat_image_colorize_field(fimg_rgba, fimg_mono, msk, zmax, style);

    /* Overlay electrode positions: */
    int32_t icop = 3; /* Opacity channel. */
    float_image_overlay(fimg_rgba, elp, icop, 0, 0);
  }

nan_options_t *nan_parse_options(int32_t argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    /* Allocate the command line argument record: */
    nan_options_t *o = notnull(malloc(sizeof(nan_options_t)), "no mem");

    /* Parse keyword parameters: */

    argparser_get_keyword(pp, "-capType");
    o->capType = argparser_get_next_non_keyword(pp);

    argparser_get_keyword(pp, "-maskDir");
    o->maskDir = argparser_get_next_non_keyword(pp);

    argparser_get_keyword(pp, "-basisDir");
    o->basisDir = argparser_get_next_non_keyword(pp);
      
    argparser_get_keyword(pp, "-inDir");
    o->inDir = argparser_get_next_non_keyword(pp);
    
    if (argparser_keyword_present(pp, "-skip"))
      { o->skip = (int32_t)argparser_get_next_int(pp, 0, INT32_MAX); }
    else
      { o->skip = 0; }
      
    if (argparser_keyword_present(pp, "-read"))
      { o->read = (int32_t)argparser_get_next_int(pp, INT32_MIN, INT32_MAX);
        if (o->read <= 0) { o->read = INT32_MIN; }
      }
    else
      { o->read = INT32_MIN; }
      
    if (argparser_keyword_present(pp, "-step"))
      { o->step = (int32_t)argparser_get_next_int(pp, 1, INT32_MAX); }
    else
      { o->step = 1; }
      
    argparser_get_keyword(pp, "-outDir");
    o->outDir = argparser_get_next_non_keyword(pp);
    
    o->setName = string_vec_new(100);
    int32_t nd = 0; /* Number of datasets specified. */
    while (argparser_keyword_present(pp, "-setName"))
      { string_vec_expand(&(o->setName), nd);
        o->setName.e[nd] = argparser_get_next_non_keyword(pp);
        nd++;
      }
    string_vec_trim(&(o->setName), nd);
    if (nd == 0) { argparser_error(pp, "must specify at least one \"-setName\""); }
    
    o->marker = nan_mkop_vec_new(100);
    int32_t nm = 0; /* Number of datasets specified. */
    while (argparser_keyword_present(pp, "-marker"))
      { nan_mkop_vec_expand(&(o->marker), nm);
        o->marker.e[nm].ic = (int32_t)argparser_get_next_int(pp, INT32_MIN, INT32_MAX);
        o->marker.e[nm].vref = argparser_get_next_double(pp, -1.0e6, 1.0e6);
        o->marker.e[nm].color.c[0] = (float)argparser_get_next_double(pp, 0.0, 1.0);
        o->marker.e[nm].color.c[1] = (float)argparser_get_next_double(pp, 0.0, 1.0);
        o->marker.e[nm].color.c[1] = (float)argparser_get_next_double(pp, 0.0, 1.0);
        nm++;
      }
    nan_mkop_vec_trim(&(o->marker), nm);

    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }
  
vec_typeimpl(nan_mkop_vec_t,nan_mkop_vec,nan_mkop_t);
