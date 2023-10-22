#define PROG_NAME "nmeeg_make_imgbasis"
#define PROG_DESC "Performs principal component analysis on an EEG dataset."
#define PROG_VERS "2013-06-06"

#define nmeeg_make_imgbasis_C_COPYRIGHT \
  "Copyright © 2013 by the State University of Campinas (UNICAMP)"
/* Last edited on 2023-10-21 21:56:09 by stolfi */
    
#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -type {BASTYPE} \\\n" \
  "    -size {NX} {NY} \\\n" \
  "    -electrodes {NE} \\\n" \
  "    -outPrefix {OUTPREFIX} \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program generates a set of {NE} images that can be combined to interpolate" \
  " a set of electrode potentials over the scalp.\n" \
  "\n" \
  "  The program currently accepts {NE=20} or {NE=128} only.  In both cases the" \
  " electrodes are assumed to be in fixed nominal position as defined" \
  " by {neuromat_eeg_geom_get_20_schematic_2D_points} and" \
  " {neuromat_eeg_geom_get_128_schematic_2D_points}.\n" \
  "\n" \
  "  All output files are sigle-channel float-valued images" \
  " as written by {float_image_write}, with the same" \
  " size, namely {NX} columns and {NY} rows of pixels.\n" \
  "\n" \
  "  The basis images will be called \"{OUTPREFIX}_b{BASTYPE}_{CHNAME[k]}.fni\" where {CHNAME[0..NE-1} are" \
  " the channel names as defined by {neuromat_eeg_get_channel_names}.  Their sample values should span" \
  " the range {[0_1]} possibly with some over- or undershoot.\n" \
  "\n" \
  "  The programs also outputs an mask image \"{OUTPREFIX}_msk.fni\" that is 1 where the" \
  " interpolation is defined and 0 where it is not, with antialiased boundary." \
  " It also outputs an antialiased overlay image \"{OUTPREFIX}_elp.fni\" that" \
  " is 1 at the electrode positions and 0 away from them, with antialiased boundaries.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -type {BASTYPE}\n" \
  "    This mandatory arguments selects the type of basis.  Currently an integer {0,1,2}.\n" \
  "\n" \
  "  -size {NX} {NY}\n" \
  "    This mandatory argument specifies the width {NX} and height {NY} of all images, in pixels.\n" \
  "\n" \
  "  -electrodes {NE}\n" \
  "    This mandatory argument specifies the number of electrodes to assume.\n" \
  "\n" \
  "  -outPrefix {OUTPREFIX}\n" \
  "    This mandatory argument specifies the common prefix of all output" \
  " images.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  neuromat_animate(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2013-06-01 as part of {nmeeg_animate.c} by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2013-12-06 J. Stolfi: split off from {nmeeg_animate.c}, substantial rewrite.\n" \
  "  2013-12-06 J. Stolfi: uses {argparser.h}.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " nmeeg_make_imgbasis_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define stringify(x) strngf(x)
#define strngf(x) #x

/* Last edited on 2013-11-30 06:15:07 by stolfilocal */
/* Converts an EEG dataset into an animation. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <r2.h>
#include <float_image.h>
#include <float_image_paint.h>
#include <affirm.h>
#include <jsfile.h>
#include <argparser.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_geom.h>
#include <neuromat_eeg_image_basis.h>

typedef struct nmib_options_t
  { int type;    /* Type of basis. */
    int size_NX; /* Width. */
    int size_NY; /* Height. */
    int electrodes;      /* Electrode count. */
    char *outPrefix;  /* Prefix for all output files. */
  } nmib_options_t;
  /* Arguments from command line. */

void nmib_image_write(char *prefix, int type, char *tag, float_image_t *fim);
  /* Writes image {fim} to file "{prefix}_b{type}_{tag}.fni", where {type}
    is typeset with two digits, zero-padded. However, omits
    the "-b{type}" part if {type} is negative. */
    
nmib_options_t *nmib_parse_options(int argc, char **argv);
  /* Parses the command line arguments. */
    
void nmib_choose_scalp_ellipse(int NX, int NY, r2_t *ctrP, r2_t *radP);
  /* Chooses the size and center of schematic head ellipse. */

int main(int argc, char **argv)
  {
    nmib_options_t *o = nmib_parse_options(argc, argv);
    
    fprintf(stderr, "output files are named %s_b%02d_*.fni\n", o->outPrefix, o->type);
    
    int NX = o->size_NX;
    int NY = o->size_NY;
    int ne = o->electrodes;
    fprintf(stderr, "assuming %d electrodes in standard positions\n", ne);
    
    /* Get the schematic electrode positions in unit disk: */
    r2_t *pos2D = neuromat_eeg_geom_get_schematic_2D_points(ne);

    /* Get the electrode names: */
    char **chname = NULL;
    int nc;
    neuromat_eeg_get_channel_names(ne, 0, NULL, &nc, &chname);
    assert(nc == ne+1);

    /* Geometry of nominal scalp projection ellipse: */
    r2_t ctr; /* Center of ellipse. */
    r2_t rad; /* Semidiameters in each axis. */
    nmib_choose_scalp_ellipse(NX, NY, &ctr, &rad);

    /* Make the head mask: */
    fprintf(stderr, "building the head outline mask\n");
    float_image_t *msk = neuromat_eeg_image_schematic_head_mask(NX, NY, &ctr, &rad);
    nmib_image_write(o->outPrefix, -1, "msk", msk);

    /* Make the channel basis images: */
    fprintf(stderr, "creating the basis images\n");
    float_image_t **bas = neuromat_eeg_image_basis_make(o->type, msk, ne, pos2D, &ctr, &rad);
    
    /* Write the channel basis images: */
    int ie;
    for (ie = 0; ie < ne; ie++) 
      { nmib_image_write(o->outPrefix, o->type, chname[ie], bas[ie]); }

    /* Map the electrode positions to image coordinates: */
    r2_t *posImg = notnull(malloc(ne*sizeof(r2_t)), "no mem");
    neuromat_eeg_geom_map_many_disk_to_ellipse(ne, pos2D, &ctr, &rad, posImg);
    
    /* Draw electrodes: */
    float_image_t *elp = float_image_new(1, NX, NY);
    double drad = 1.50;       /* Dot radius. */
    float bck = 0.0;
    float wht = 1.0;
    neuromat_eeg_image_draw_electrodes(elp, ne, posImg, drad, -1, &bck, &wht);
    nmib_image_write(o->outPrefix, -1, "elp", elp);
        
    float_image_free(msk);
    float_image_free(elp);
    { int ie; for (ie = 0; ie < ne; ie++) { float_image_free(bas[ie]); } }
    free(bas);
    
    return 0;
  }

void nmib_image_write(char *prefix, int type, char *tag, float_image_t *fim)
  { 
    char *fname = NULL;
    if (type < 0)
      { asprintf(&fname, "%s_%s.fni", prefix, tag); }
    else
      { asprintf(&fname, "%s_b%02d_%s.fni", prefix, type, tag); }
    FILE *wr = open_write(fname, TRUE);
    float_image_write(wr, fim);
    fclose(wr);
    free(fname);
  }
  
void nmib_choose_scalp_ellipse(int NX, int NY, r2_t *ctrP, r2_t *radP)
  {
    double aspect = 1.00; /* Width over height. */
    double size = 0.95*fmin(NX/aspect, NY); /* Height of ellipse. */
    r2_t ctr = (r2_t){{ 0.50*NX, 0.50*NY }}; (*ctrP) = ctr;
    r2_t rad = (r2_t){{ aspect*size/2, size/2 }}; (*radP) = rad;
  }

nmib_options_t *nmib_parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    /* Allocate the command line argument record: */
    nmib_options_t *o = notnull(malloc(sizeof(nmib_options_t)), "no mem");

    /* Parse keyword parameters: */

    argparser_get_keyword(pp, "-type");
    o->type = (int)argparser_get_next_int(pp, 0, INT_MAX);

    argparser_get_keyword(pp, "-size");
    o->size_NX = (int)argparser_get_next_int(pp, 1, INT_MAX);
    o->size_NY = (int)argparser_get_next_int(pp, 1, INT_MAX);

    argparser_get_keyword(pp, "-electrodes");
    o->electrodes = (int)argparser_get_next_int(pp, 1, INT_MAX);
    
    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next_non_keyword(pp);

    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }
