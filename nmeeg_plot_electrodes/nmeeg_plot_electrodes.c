#define PROG_NAME "nmeeg_plot_electrodes"
#define PROG_DESC "Plots the electrode positions for the NeuroMat project"
#define PROG_VERS "1.0"

#define nmeeg_plot_electrodes_C_COPYRIGHT \
  "Copyright © 2013 by the State University of Campinas (UNICAMP)"

/* Last edited on 2023-10-22 19:38:10 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -capType {CAP_TYPE} \\\n" \
  "    -radius {RADIUS} \\\n" \
  "    -outPrefix {OUT_PREFIX} \\\n" \
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
  "  -capType {CAP_TYPE}\n" \
  "    This mandatory argument specifies the cap type, which determines" \
  " the count, names, and positions of the electrodes.  The possible" \
  " values of {CAP_TYPE} are defined" \
  " by {neuromat_eeg_get_channel_names} in {neuromat_eeg.h}, and" \
  " include \"R20\", \"R128\", \"R129\", and \"FN3\".\n" \
  "\n" \
  "  -radius {RADIUS} \\\n" \
  "    This mandatory argument specifies the radius (in mm) of the circle" \
  " that represents the subject's head.\n" \
  "\n" \
  "  -outPrefix {OUT_PREFIX} \\\n" \
  "    This mandatory argument specifies the common prefix for all file names. The files" \
  " will be called \"{OUT_PREFIX}-{NUM}.eps\" where {NUM} is a sequental number.\n" \
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
  "  2023-10-22 J. Stolfi: replaced {pswr.h} by {epswr.h}.\n" \
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

#include <epswr.h>
#include <jsstring.h>
#include <affirm.h>
#include <argparser.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_geom.h>

/* General plot options: */
typedef struct npe_options_t
  { char *capType;    /* Type of EEG cap, like "R20", "R128" etc. */
    char *outPrefix;    /* Prefix for output file names. */
    double radius;    /* Nominal head circle radius (in mm). */
  } npe_options_t;

/* INTERNAL PROTOTYPES */

int32_t main(int32_t argc, char **argv);

npe_options_t *npe_parse_options(int32_t argc, char **argv);

void npe_plot_electrodes(epswr_figure_t *eps, npe_options_t *o);
  /* Plots the desired diagram to the Postscript stream {ps}. */

epswr_figure_t *npe_start_figure(npe_options_t *o, char *figTag, double radius);
  /* Creates a new EPS figure writer object {eps} whose filename is
    "{o->outPrefix}-{figTag}.eps". The nominal client plot window will be
    {[-radius _ +radius]^2} plus a small margin. */

void npe_finish_figure(epswr_figure_t *eps, npe_options_t *o);
  /* Finalizes the current figure written to file {eps}. */

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {   
    npe_options_t *o = npe_parse_options(argc, argv);
    
    epswr_figure_t *eps = npe_start_figure(o, "els", o->radius);
    npe_plot_electrodes(eps, o);
    npe_finish_figure(eps, o);
    return(0);
  }

void npe_plot_electrodes(epswr_figure_t *eps, npe_options_t *o)
  {
    int32_t ne = -1;
    char **chname = NULL;
    r2_t *pos2D = NULL;
    neuromat_eeg_geom_get_schematic_2D_points(o->capType, &ne, &chname, &pos2D);
    
    epswr_set_pen(eps, 0,0,0, 0.25, 0,0);
    epswr_circle(eps, 0.0, 0.0, o->radius, FALSE, TRUE);
    
    epswr_set_fill_color(eps, 0.800, 1.000, 0.500);
    epswr_set_pen(eps, 0,0,0, 0.10, 0,0);
    epswr_set_label_font(eps, "CourierBold", 8);
    
    for (int32_t i = 0; i < ne; i++)
      { /* Grab the electrode {i}: */
        r2_t *pi = &(pos2D[i]);
        double xi = o->radius * pi->c[0];
        double yi = o->radius * pi->c[1];
        char *cni = chname[i];
        /* Suppress the 'C' in the electrode names of the 128-electrode setup: */
        if ((ne == 128) && ((*cni) == 'C')) { cni++; }
        /* Plot the electrode: */
        epswr_dot(eps, xi, yi, 3.75, TRUE, TRUE);
        epswr_label(eps, cni, "R", xi, yi, 0.0, TRUE, 0.5,0.5, TRUE, TRUE);
      }
  }

epswr_figure_t *npe_start_figure(npe_options_t *o, char *figTag, double radius)
  { 
    double mm = epswr_pt_per_mm;
    double mrg = 5.0; /* Margin (pt) */
    
    double xymin = -1.01*radius;
    double xymax = +1.01*radius;

    double hPlotSize = (xymax - xymin) * mm;
    double vPlotSize = (xymax - xymin) * mm;
    
    epswr_figure_t *eps = epswr_new_named_figure
      ( NULL, o->outPrefix, figTag, -1, NULL, 
        hPlotSize, vPlotSize,
        mrg, mrg, mrg, mrg,
        TRUE
      );
      
    epswr_set_client_window(eps, xymin, xymax, xymin, xymax);
    return eps;
  }
  
void npe_finish_figure(epswr_figure_t *eps, npe_options_t *o)
  { epswr_end_figure(eps); }
  
#define npe_MIN_RADIUS (10.0)
#define npe_MAX_RADIUS (200.0)
#define npe_MAX_NELEC (128)

npe_options_t *npe_parse_options(int32_t argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    npe_options_t *o = notnull(malloc(sizeof(npe_options_t)), "no mem"); 
    
    /* Parse keyword parameters: */
 
    argparser_get_keyword(pp, "-capType");
    o->capType = argparser_get_next_non_keyword(pp);
   
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

