#define PROG_NAME "nmeeg_make_imgbasis"
#define PROG_DESC "Generates an interpolation/approximation basis for a set of electrodes."
#define PROG_VERS "2013-06-06"

#define nmeeg_make_imgbasis_C_COPYRIGHT \
  "Copyright © 2013 by the State University of Campinas (UNICAMP)"
/* Last edited on 2023-10-22 09:47:08 by stolfi */
    
#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -capType {CAP_TYPE} \\\n" \
  "    -basisType { shepard | gauss | mexican | voronoi } \\\n" \
  "    -interp { F | T } \\\n" \
  "    -norm { F | T } \\\n" \
  "    -scalpSize {SX} {SY} {SZ} \\\n" \
  "    -imageSize {NX} {NY} \\\n" \
  "    -aspect {ASPECT} \\\n" \
  "    -extent {EXTENT} \\\n" \
  "    [ -subsample {MSUB} ] \\\n" \
  "    -outDir {OUTDIR} \\\n" \
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
  "  The program currently assumes that the electrodes are a subset of the standard 20-electrode" \
  " or 129-electrode layouts used at INDC/UFRJ.  In each case the electrodes are assumed" \
  " to located at fixed /schematic 2D positions/ as defined" \
  " by {neuromat_eeg_geom_get_schematic_2D_points}.  They are therfore also in" \
  " fixed /schematic 3D positions/ as defined by {neuromat_eeg_geom_3D_from_2D}.\n" \
  "\n" \
  "  All basis image fields are sigle-channel float-valued images" \
  " as written by {float_image_write}, with the same" \
  " size, namely {NX} columns and {NY} rows of pixels.\n" \
  "\n" \
  "  The basis images will be called \"{OUTDIR}/{IMGDIR}/{BTYPE}_i{INTERP}_n{NORM}/{CHNAME[k]}.fni\" where" \
  " {CHNAME[0..NE-1} are the channel names as defined by {neuromat_eeg_get_channel_names}; and {IMGDIR}" \
  " is \"{NX}_{NY}\", the numbers being padded with zeros to 4, 4, and 3 digits. The" \
  " sample values generally should span the range {[0_1]} possibly with some over- or undershoot.\n" \
  "\n" \
  "  The programs also outputs an mask image \"{OUTDIR}/{IMGDIR}/msk.fni\" that is 1 where the" \
  " interpolation is defined and 0 where it is not, with antialiased boundary.\n" \
  "\n" \
  " Finally the program also outputs an antialiased overlay RGBA image \"{OUTDIR}/{IMGDIR}/elp.fni\" that" \
  " has the electrode postions marked as dots, with antialiased boundaries, over a" \
  " transparent background.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -capType {CAP_TYPE}\n" \
  "    This mandatory argument specifies the cap type, which determines" \
  " the count, names, and positions of the electrodes.  The possible" \
  " values of {CAP_TYPE} are defined" \
  " by {neuromat_eeg_get_channel_names} in {neuromat_eeg.h}, and" \
  " include \"R20\", \"R128\", \"R129\", and \"FN3\".\n" \
  "\n" \
  "  -basisType { shepard | gauss | mexican | voronoi }\n" \
  "    This mandatory argument selects the type of basis.\n" \
  "\n" \
  "  -interp { F | T } \n" \
  "     This mandatory argument specifies whether the basis should be strictly" \
  " interpolating.  Namely, whether the interpolated value at each electrode position" \
  " should be the measured value at that electrode. Note that basis" \
  " types \"shepard\" and \"vorono\" are always interpolating.\n" \
  "\n" \
  "  -norm { F | T } \n" \
  "     This mandatory argument specifies whether the values of the basis elements" \
  " at each point {p} should be scaled so that their sum is 1.  This property is" \
  " required for the basis toreproduce constant functions; that is, so that the" \
  " interpolated vaile at every point be the same value {V} when all the measured" \
  " electrode potentials are equa to {V}.  Note that basis" \
  " types \"shepard\" and \"vorono\" always have this property." \
  "\n" \
  "  -scalpSize {SX} {SY} {SZ}\n" \
  "    This mandatory argument specifies the width {SX}, the length {SY}, and the" \
  " height {SZ} (in millimeters) of the idealized scalp.  It affects the assumed" \
  " electrode positions in 3D, and hence the distances between electrode positions" \
  " and  other points of the scalp.  The ratio of {SX} to {SY} also determines the" \
  " aspect ratio of the 2D scalp on the images.\n" \
  "\n" \
  "  -imageSize {NX} {NY}\n" \
  "    This mandatory argument specifies the width {NX} and height {NY} of all images, in pixels.\n" \
  "\n" \
  "  -extent {EXTENT}\n" \
  "    This mandatory argument specifies the maximum extent of the scalp ellipse relative" \
  " to the image size, along any axis.  Thus \"-extent 1.00\" will result in the largest scalp" \
  " ellipse of the specified aspect ratio that fits in the image.\n" \
  "\n" \
  "  -subsample {MSUB}\n" \
  "    This optional argument specifies the degre of sub-sampling to use.  The" \
  " value of each pixel in each basis image will be the average of the corresponding" \
  " basis function values at a grid of {MSUB} by {MSUB} points within the pixel.  If" \
  " omitted, the program assumes \"-subsample 1\". Note that subsampling" \
  " has a discernible effect only for basis functions that are not" \
  " continuous, like \"voronoi\".\n" \
  "\n" \
  "  -outDir {OUTDIR}\n" \
  "    This mandatory argument specifies the diretory where all output files will be placed." \
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
  "  2017-09-27 J. Stolfi: added \"-aspect\" and \"-extent\" options.\n" \
  "  2021-08-22 J. Stolfi: added \"-interp\" and \"-norm\" options. Adapted to new library.\n" \
  "  2021-08-22 J. Stolfi: added \"-scalpSize\" and \"-subsample\" options.\n" \
  "  2021-08-28 J. Stolfi: replaced \"-outPrefix\" by \"-outDir\".\n" \
  "  2021-08-31 J. Stolfi: changed \"-electrodes\" from count to list of names.\n" \
  "  2023-10-22 J. Stolfi: replaced \"-electrodes\" by \"-capType\".\n" \
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
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <r2.h>
#include <vec.h>
#include <float_image.h>
#include <float_image_paint.h>
#include <affirm.h>
#include <jsfile.h>
#include <argparser.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_geom.h>
#include <neuromat_eeg_image_basis.h>

typedef enum {
    nmib_basis_type_SHEPARD,
    nmib_basis_type_GAUSS,
    nmib_basis_type_MEXICAN,
    nmib_basis_type_VORONOI
  } nmib_basis_type_t;
  /* The avaliable basis types. */
  
#define nmib_basis_type_MAX nmib_basis_type_VORONOI
  /* Max {nmib_basis_type_t} value. */
  
typedef struct nmib_options_t
  { char *capType;              /* Type of EEG cap, like "R20", "R128" etc. */
    nmib_basis_type_t basisType;     /* Type of basis. */
    bool_t interp;              /* Whether the basis is interpolatory. */
    bool_t norm;                /* Whether the basis values should add to 1 at each point. */
    r3_t scalpSize;             /* Diameter of the idealized 3D scalp ellipsoid along each axis. */
    int32_t imageSize_NX;       /* Width. */
    int32_t imageSize_NY;       /* Height. */
    double extent;              /* Ratio of idealized 2D scalp ellipse size to the max size that fits in image. */
    int32_t subsample;          /* Subsampling order in each pixel. */
    char *outDir;               /* Directory for all output files. */
  } nmib_options_t;
  /* Arguments from command line. */
    
float_image_t **nmib_make_basis
  ( nmib_basis_type_t btype,        /* Type of basis. */
    bool_t interp,                  /* Whether the basis is interpolatory. */
    bool_t norm,                    /* Whether the basis values should add to 1 at each point. */
    int32_t ne,                     /* Electrode count. */
    r2_t pos2D[],                   /* Schematic 2D electrode positions. */    
    r3_t *srad,                     /* Radius of the idealized 3D scalp ellipsoid along each axis. */
    float_image_t *msk,             /* Scalp domain mask. */
    r2_t *ictr,                     /* Center of the idealized 2D skull ellipse on image. */
    r2_t *irad,                     /* Radii of the idealized 2D skull ellipse in each axis. */
    int32_t msub                    /* Order of subsampling of the basis functions in each pixel. */
  );
  /* Returns ar array of {ne} images {img[0..ne-1]} where {img{k]} shows
    the influence of the potential measured at electrode {k} on the
    interpolated/approximated potential at each point if the scalp.
    
    The image points are mapped to /schematic 2D coordinates/ by the
    affine map that takes the ellipse with center {ictr} and radii
    {irad} to the canonical unit disk. They are then mapped to points of
    the idealized 3D skull, the canonical ellipsoid with radii {srad}.
    
    A pixel of a basis image is computed only is the corresponding pixel
    of the mask {msk} is nonzero. Each pixel value will be the averahe
    of the corresponding basis function value at the points on the
    idealized skull that correspond to a grid of {msub} times {msub}
    points inside the pixel. */
    
char *const nmib_basis_name(nmib_basis_type_t btype);
  /* Returns the name of the basis with type {btype}. */

void nmib_image_write(char *imgDir, int32_t btype, bool_t interp, bool_t norm, char *tag, float_image_t *fim);
  /* Writes image {fim} to a ".fni" image file.
    
    If {btype} is non-negative, it must be a valid {nmib_basis_type_t}
    value; the file is then named
    "{prefix}_{bname}_i{interp}_n{norm}_{tag}.fni" where {bname} is the
    basis name. The booleans {interp} and {norm} are converted to 0 and
    1.
  
    If {btype} is negative, the file will be named "{prefix}_{tag}.fni". The parameters 
    {interp} and {norm} are ignored. */
    
nmib_options_t *nmib_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments. */
    
void nmib_choose_scalp_ellipse(int32_t NX, int32_t NY, double aspect, double extent, r2_t *ctrP, r2_t *radP);
  /* Chooses the size and center of schematic scalp ellipse. 
    
    The ellipse will be centered in the rectangle {[0_NX]×[0_NY]}.  
    The {aspect} parameter is the desired ratio of width to heightof the ellipse.  
    The size of the ellipse will be {extent} times the size of the
    largest ellipse of the given aspect that fits in the image.  */

int32_t main(int32_t argc, char **argv)
  {
    nmib_options_t *o = nmib_parse_options(argc, argv);
    
    int32_t NX = o->imageSize_NX;
    int32_t NY = o->imageSize_NY;

    /* Get the schematic electrode positions in unit disk: */
    int32_t ne = -1;  /* Number of electrodes. */
    char **chname = NULL;  /* Channel names. */
    neuromat_eeg_get_channel_names(o->capType, 0, NULL, &ne, &chname);
    r2_t *pos2D = neuromat_eeg_geom_get_schematic_2D_points_by_name(o->capType, ne, chname);

    /* Radii of the idealized 3D scalp ellipsoid: */
    r3_t srad = (r3_t){{ o->scalpSize.c[0]/2,  o->scalpSize.c[1]/2,  o->scalpSize.c[2]/2 }};

    /* Geometry of idealized 2D scalp projection ellipse on image: */
    r2_t ictr; /* Center of ellipse. */
    r2_t irad; /* Semidiameters in each axis. */
    double aspect = srad.c[0] / srad.c[1];
    nmib_choose_scalp_ellipse(NX, NY, aspect, o->extent, &ictr, &irad);
    fprintf(stderr, "upper scalp in projection:\n");
    fprintf(stderr, "   ctr = ( %.1f %.1f )\n", ictr.c[0], ictr.c[1]);
    fprintf(stderr, "   rad = ( %.1f %.1f )\n", irad.c[0], irad.c[1]);

    /* Subdiretory for basis files: */
    char *imgDir = NULL;
    asprintf(&imgDir, "%s/%04d_%04d", o->outDir, NX, NY);
    
    /* Make the head mask: */
    fprintf(stderr, "building the head outline mask\n");
    float_image_t *msk = neuromat_eeg_image_make_idealized_scalp_mask(NX, NY, &ictr, &irad);
    nmib_image_write(imgDir, -1, FALSE, FALSE, "msk", msk);

    /* Make the channel basis images: */
    fprintf(stderr, "creating the basis images\n");
    int32_t msub = o->subsample;
    float_image_t **bas = nmib_make_basis(o->basisType, o->interp, o->norm, ne, pos2D, &srad, msk, &ictr, &irad, msub);
    
    /* Write the channel basis images: */
    int32_t ie;
    for (ie = 0; ie < ne; ie++) 
      { nmib_image_write(imgDir, o->basisType, o->interp, o->norm, chname[ie], bas[ie]); }

    /* Map the electrode positions to image coordinates: */
    r2_t *posImg = notnull(malloc(ne*sizeof(r2_t)), "no mem");
    neuromat_eeg_geom_map_many_disk_to_ellipse(ne, pos2D, &ictr, &irad, posImg);
    
    /* Draw electrodes: */
    double drad = 3.25;   /* Dot radius (px). */
    double hwd = 0.90;    /* Half-width of dot outline (px). */
    frgb_t yel = (frgb_t){{ 1.000f, 0.950f, 0.000f }};
    frgb_t blu = (frgb_t){{ 0.000f, 0.050f, 1.000f }};
    float_image_t *elp = neuromat_eeg_image_electrodes_overlay(ne, posImg, drad, hwd, -1, &yel, &blu, NX, NY);
    nmib_image_write(imgDir, -1, FALSE, FALSE, "elp", elp);
        
    float_image_free(msk);
    float_image_free(elp);
    { int32_t ie; for (ie = 0; ie < ne; ie++) { float_image_free(bas[ie]); } }
    free(bas);
    
    return 0;
  }

float_image_t **nmib_make_basis
  ( nmib_basis_type_t btype,   /* Type of basis. */
    bool_t interp,      /* Whether the basis is interpolatory. */
    bool_t norm,        /* Whether the basis values should add to 1 at each point. */
    int32_t ne,         /* Electrode count. */
    r2_t pos2D[],       /* Schematic 2D electrode positions. */    
    r3_t *srad,         /* Radius of the idealized 3D scalp ellipsoid along each axis. */
    float_image_t *msk, /* Scalp domain mask. */
    r2_t *ictr,         /* Center of the idealized 2D skull ellipse on image. */
    r2_t *irad,         /* Radii of the idealized 2D skull ellipse in each axis. */
    int32_t msub        /* Order of subsampling of the basis functions in each pixel. */
  )
  {
    /* Map schematic 2D electrode positions to idealized 3D scalp: */
    r3_t pos3D[ne];
    for (int32_t ie = 0; ie < ne; ie++) 
      { pos3D[ie] = neuromat_eeg_geom_3D_from_2D(&(pos2D[ie]), NULL, srad); }
    
    auto void mother(int32_t ne, double bval[], r3_t *p3D);
      /* Procedure that fills {bval[0..ne-1]} with the values of the mother basis functions at point {p}. */
  
    /* Compute the idealized 3D influence radius {erad[ie]} for each electrode, if needed: */
    bool_t mother_needs_erad = 
      (btype == nmib_basis_type_SHEPARD) ||
      (btype == nmib_basis_type_GAUSS) ||
      (btype == nmib_basis_type_MEXICAN);
    double *erad = (mother_needs_erad ? neuromat_eeg_func_basis_nearest_dists(ne, pos3D) : NULL);

    /* Compute the Lagrange interpolation matrix {L} that converts the mother basis to interpolating basis: */
    bool_t mother_is_interp =
      (btype == nmib_basis_type_SHEPARD) || 
      (btype == nmib_basis_type_VORONOI);
    double *L = (interp && (! mother_is_interp) ? neuromat_eeg_func_basis_lagrangian_matrix(ne, mother, pos3D) : NULL);
    
    auto void eval(int32_t ne, double bval[], r3_t *p3D);
      /* Procedure that fills {bval[0..ne-1]} with the values of the final basis functions at point {p}. */
    
    float_image_t **img = neuromat_eeg_image_basis_make(ne, eval, msk, msub, srad, ictr, irad);
    
    return img;
    
    /* Internal procedures: */
    
    auto void mother(int32_t ne, double bval[], r3_t *p3D)
      {
        for (int32_t ie = 0; ie < ne; ie++) 
          { 
            switch(btype)
              { 
                case nmib_basis_type_SHEPARD:
                  { double rho = erad[ie];
                    double sigma = 4*erad[ie];
                    bval[ie] = neuromat_eeg_func_basis_shepard_weight(p3D, &(pos3D[ie]), rho, sigma, 2); }
                  break;
                case nmib_basis_type_GAUSS:
                  { double sigma = erad[ie];
                    bval[ie] = neuromat_eeg_func_basis_gauss_bell(p3D, &(pos3D[ie]), sigma);  }
                  break;
                case nmib_basis_type_MEXICAN:
                  { double sigma = 4*erad[ie];
                    double tau = 2*erad[ie];
                    bval[ie] = neuromat_eeg_func_basis_mexican_hat(p3D, &(pos3D[ie]), sigma, tau); }
                  break;
                case nmib_basis_type_VORONOI:
                  { bval[ie] = neuromat_eeg_func_basis_voronoi_ind(p3D, ie, ne, pos3D); }
                  break;
                default:
                  assert(FALSE);
              }
          }
      }
  
    void eval(int32_t ne1, double bval[], r3_t *p3D)
      { assert(ne1 == ne);
        neuromat_eeg_func_basis_eval(ne, bval, p3D, mother, L, norm);
      }
  }
 
char *const nmib_basis_name(nmib_basis_type_t btype)
  { switch(btype)
      { 
        case nmib_basis_type_SHEPARD:
          return "shepard";
        case nmib_basis_type_GAUSS:
          return "gauss";
        case nmib_basis_type_MEXICAN:
          return "mexican";
        case nmib_basis_type_VORONOI:
          return "voronoi";
        default:
          assert(FALSE);
      }
  }

void nmib_image_write(char *imgDir, int32_t btype, bool_t interp, bool_t norm, char *tag, float_image_t *fim)
  { 
    char *fname = NULL;
    if (btype < 0)
      { asprintf(&fname, "%s/%s.fni", imgDir, tag); }
    else
      { char *const bname = nmib_basis_name((nmib_basis_type_t)btype);
        asprintf(&fname, "%s/%s_i%d_n%d/%s.fni", imgDir, bname, interp, norm, tag);
      }
    FILE *wr = open_write(fname, TRUE);
    float_image_write(wr, fim);
    fclose(wr);
    free(fname);
  }
  
void nmib_choose_scalp_ellipse(int32_t NX, int32_t NY, double aspect, double extent, r2_t *ctrP, r2_t *radP)
  {
    double maxrad = fmin(0.50*NX/aspect, 0.50*NY);
    r2_t ictr = (r2_t){{ 0.50*NX, 0.50*NY }};
    r2_t irad = (r2_t){{ extent*aspect*maxrad, extent*maxrad }};
    (*ctrP) = ictr;
    (*radP) = irad;
  }

nmib_options_t *nmib_parse_options(int32_t argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    /* Allocate the command line argument record: */
    nmib_options_t *o = notnull(malloc(sizeof(nmib_options_t)), "no mem");

    /* Parse keyword parameters: */

    argparser_get_keyword(pp, "-capType");
    o->capType = argparser_get_next_non_keyword(pp);

    argparser_get_keyword(pp, "-basisType");
    char *bname = argparser_get_next_non_keyword(pp);
    o->basisType = -1;
    for (int32_t btype = 0; btype <= nmib_basis_type_MAX; btype++)
      { if (strcmp(bname, nmib_basis_name(btype)) == 0)
          { o->basisType = (nmib_basis_type_t)btype; break; }
      }
    if (o->basisType < 0) { argparser_error(pp, "invalid basis type"); }

    argparser_get_keyword(pp, "-interp");
    o->interp = argparser_get_next_bool(pp);

    argparser_get_keyword(pp, "-norm");
    o->norm = argparser_get_next_bool(pp);

    argparser_get_keyword(pp, "-scalpSize");
    for (int32_t i = 0; i < 3; i++) 
      { o->scalpSize.c[i] = argparser_get_next_double(pp, 1.0, 300.0); }

    argparser_get_keyword(pp, "-imageSize");
    o->imageSize_NX = (int32_t)argparser_get_next_int(pp, 1, INT32_MAX);
    o->imageSize_NY = (int32_t)argparser_get_next_int(pp, 1, INT32_MAX);

    argparser_get_keyword(pp, "-extent");
    o->extent = argparser_get_next_double(pp, 0.01, 1.00);
    
    if (argparser_keyword_present(pp, "-subsample"))
      { o->subsample = (int32_t)argparser_get_next_int(pp, 1, 5); }
    else
      { o->subsample = 1; }
    
    argparser_get_keyword(pp, "-outDir");
    o->outDir = argparser_get_next_non_keyword(pp);

    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }
