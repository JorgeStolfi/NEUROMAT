#define PROG_NAME "nmeeg_plot_electrodes"
#define PROG_DESC "Plots the electrode positions for the NeuroMat project"
#define PROG_VERS "1.0"

#define nmeeg_plot_electrodes_C_COPYRIGHT \
  "Copyright © 2013 by the State University of Campinas (UNICAMP)"

/* Last edited on 2023-10-22 11:05:54 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -capType {CAP_TYPE} \\\n" \
  "    -radius {RADIUS} \\\n" \
  "    -outPrefix {OUTNAME} \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program writes an encapsulated Postscript plot of the nominal electrode positions.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -radius {RADIUS} \\\n" \
  "    This mandatory argument specifies the radius (in mm) of the circle" \
  " that represents the subject's head.\n" \
  "\n" \
  "  -nElec {NELEC} \\\n" \
  "    This mandatory argument specifies the number of electrodes.  Currently" \
  " only {NELEC=20} and {NELEC=128} are allowed.\n" \
  "\n" \
  "  -outPrefix {OTNAME} \\\n" \
  "    This mandatory argument specifies the common prefix for all file names. The files" \
  " will be called \"{OUTNAME}-{NUM}.eps\" where {NUM} is a sequental number.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  what is in front of you.\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2013-09-30 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2023-10-22 J. Stolfi: replaced \"-nElec\" by \"-capType\".\n" \
  "  2023-10-22 J. Stolfi: replaced \"-outName\" by \"-outPrefix\".\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " nmeeg_plot_electrodes_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <values.h>

#include <pswr.h>
#include <jsstring.h>
#include <affirm.h>
#include <argparser.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_geom.h>

/* General plot options: */
typedef struct npe_options_t
  { int nElec;        /* Number of bona-fide electrodes (excl. markers, ref volt). */
    char *outPrefix;    /* Prefix for output file names. */
    double radius;    /* Nominal head circle radius (in mm). */
  } npe_options_t;

/* INTERNAL PROTOTYPES */

int main(int argc, char **argv);

npe_options_t *npe_parse_options(int argc, char **argv);

void npe_plot_electrodes(PSStream *ps, npe_options_t *o);
  /* Plots the desired diagram to the Postscript stream {ps}. */

PSStream *npe_start_document(npe_options_t *o);
  /* IF {o->eps} is FALSE, sets {ps} to a new Postscript document
    whose filename is "{o->outPrefix}.ps".  If {o->eps} is TRUE,
    has no effect. */

void npe_start_figure(PSStream *ps, npe_options_t *o, char *figTag, double radius);
  /* The stream {ps} must be NULL; opens a new Encapsulated
    Postscript file, with name "{o->outPrefix}-{figTag}.eps", and
    initializes it.
    
    In either case the nominal client plot window is {[-radius _ +radius]^2}
    plus a small margin. */

void npe_finish_figure(PSStream *ps, npe_options_t *o);
  /* Finalizes the current figure written to file {ps}.
    If {o->eps} is true, closes the file and returns NULL.
    Otherwise simply terminates the current page. */

void npe_finish_document(PSStream *ps, npe_options_t *o);
  /* If {o->eps} is true, {ps} must be NULL and the call is 
    a no-op.  Otherwise {ps} must be a plain Postscript file;
    finalizes and closes it. */

/* IMPLEMENTATIONS */

int main(int argc, char **argv)
  {   
    npe_options_t *o = npe_parse_options(argc, argv);
    
    PSStream *ps = npe_start_document(o);
    npe_start_figure(ps, o, "els", o->radius);
    npe_plot_electrodes(ps, o);
    npe_finish_figure(ps, o);
    npe_finish_document(ps, o);
    return(0);
  }

void npe_plot_electrodes(PSStream *ps, npe_options_t *o)
  {
    int ne = o->nElec;
    int nc;
    char **chname = NULL;
    neuromat_eeg_get_channel_names(ne, 0, NULL, &nc, &chname);
    r2_t *pos = neuromat_eeg_geom_get_schematic_2D_points(ne);
    
    pswr_set_pen(ps, 0,0,0, 0.25, 0,0);
    pswr_circle(ps, 0.0, 0.0, o->radius, FALSE, TRUE);
    
    pswr_set_fill_color(ps, 0.800, 1.000, 0.500);
    pswr_set_pen(ps, 0,0,0, 0.10, 0,0);
    pswr_set_label_font(ps, "CourierBold", 8);
    
    int i;
    for (i = 0; i < ne; i++)
      { /* Grab the electrode {i}: */
        r2_t *pi = &(pos[i]);
        double xi = o->radius * pi->c[0];
        double yi = o->radius * pi->c[1];
        char *cni = chname[i];
        /* Suppress the 'C' in the electrode names of the 128-electrode setup: */
        if ((ne == 128) && ((*cni) == 'C')) { cni++; }
        /* Plot the electrode: */
        pswr_dot(ps, xi, yi, 3.75, TRUE, TRUE);
        pswr_label(ps, cni, xi, yi, 0.5, 0.5);
      }
  }

PSStream *npe_start_document(npe_options_t *o)
  { char *prefix = txtcat(o->outPrefix, "-");
    PSStream *ps = pswr_new_stream(prefix, NULL, TRUE, NULL, NULL, 640.0, 480.0);
    return ps;
  }

void npe_start_figure(PSStream *ps, npe_options_t *o, char *figTag, double radius)
  { 
    fprintf(stderr, "=== %s =========================\n", figTag);

    double margin[2]; /* X and Y extra margin (in mm). */
    margin[0] = 1.0;
    margin[1] = 1.0;
    
    double xsz = 2.01*radius;
    double ysz = 2.01*radius;
    
    interval_t bbox[2];  /* Client plot window. */
    bbox[0] = (interval_t){{ -xsz/2, +xsz/2 }};
    bbox[1] = (interval_t){{ -ysz/2, +ysz/2 }};
    
    double pt_per_mm = 72.0/25.4;

    /* Useful figure size in pt: */
    double hsz = xsz * pt_per_mm;
    double vsz = ysz * pt_per_mm;

    /* Figure margin width in pt: */
    double hmg = margin[0] * pt_per_mm;
    double vmg = margin[1] * pt_per_mm;

    /* Total figure size in pt: */
    double htot = hsz + 2*hmg;
    double vtot = vsz + 2*vmg;
    
    affirm(ps != NULL, "no output stream");
    
    /* Assumes encapsulated Postscript: */
    pswr_set_canvas_size(ps, htot, vtot);
    pswr_new_canvas(ps, figTag);
    pswr_set_canvas_layout(ps, hsz, vsz, FALSE, hmg, vmg, 0, 1, 1);
    pswr_new_picture(ps, LO(bbox[0]), HI(bbox[0]), LO(bbox[1]), HI(bbox[1]));
  }
  
void npe_finish_figure(PSStream *ps, npe_options_t *o)
  {  }

void npe_finish_document(PSStream *ps, npe_options_t *o)
  { if (ps != NULL) 
      { pswr_close_stream(ps); }
  }
  
#define npe_MIN_RADIUS (10.0)
#define npe_MAX_RADIUS (200.0)
#define npe_MAX_NELEC (128)

npe_options_t *npe_parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    npe_options_t *o = notnull(malloc(sizeof(npe_options_t)), "no mem"); 
    
    /* Parse keyword parameters: */
    
    argparser_get_keyword(pp, "-nElec");
    o->nElec = (int)argparser_get_next_double(pp, 0, npe_MAX_NELEC);
    
    argparser_get_keyword(pp, "-radius");
    o->radius = argparser_get_next_double(pp, npe_MIN_RADIUS, npe_MAX_RADIUS);
    
    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next(pp);
    
    /* Parse positional parameters: */
    argparser_skip_parsed(pp);

    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }

