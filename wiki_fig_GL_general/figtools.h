/* Tools for generating figures for IA/AA papers. */
/* Last edited on 2023-12-08 08:20:45 by stolfi */

#ifndef figtools_H
#define figtools_H

#include <stdio.h>
#include <values.h>

#include <epswr.h>
#include <frgb.h>
#include <interval.h>

#ifndef INF
#define INF INFINITY
  /* Infinity in double format. */
#endif

/* General plot options: */
typedef struct figtools_options_t
  { char *outName;    /* Prefix for output file names. */
    bool_t color;     /* TRUE uses colors, FALSE uses only grays. */
    double scale;     /* Figure scale factor (client to mm). */
    double margin[2]; /* X and Y extra margin (in mm). */
  } figtools_options_t;

epswr_figure_t *figtools_start_figure
  ( figtools_options_t *o,
    char *figTag,
    interval_t bbox[2]      /* Client plot window (client units). */
  );
  /* Opens a new Encapsulated Postscript file, with name 
    "{o->outName}-{figTag}.eps", and initializes it.
    Returns the {epswr_figure_t}.
    
    The nominal client plot window is {bbox}, plus a 
    safety margin {o->margin}. */

void figtools_finish_figure(epswr_figure_t *eps, figtools_options_t *o);
  /* Finalizes the current figure written to file {eps}.
    Closes the file and returns NULL.
    Otherwise simply terminates the current page. */

void figtools_hatch_rectangle(epswr_figure_t *eps, interval_t B[], double hstep);
  /* Draws hatch-lines across the cell {B}, spaced {hstep} in client
    units, with the current pen. */

#endif
