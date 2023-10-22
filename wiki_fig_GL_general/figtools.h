/* Tools for generating figures for IA/AA papers. */
/* Last edited on 2019-10-19 02:57:51 by jstolfi */

#ifndef figtools_H
#define figtools_H

#include <stdio.h>
#include <values.h>

#include <pswr.h>
#include <frgb.h>
#include <interval.h>

#ifndef INF
#define INF INFINITY
  /* Infinity in double format. */
#endif

#define mm (72.27/25.4)

/* General plot options: */
typedef struct PlotOptions
  { char *outName;    /* Prefix for output file names. */
    bool_t color;     /* TRUE uses colors, FALSE uses only grays. */
    double scale;     /* Figure scale factor (client to mm). */
    double margin[2]; /* X and Y extra margin (in mm). */
  } PlotOptions;

PSStream *start_figure
  ( PlotOptions *o,
    char *figTag,
    interval_t bbox[2],    /* Client plot window. */
    double scale            /* Plot scale (mm per client unit). */
  );
  /* Opens a new Encapsulated Postscript file, with name 
    "{o->outName}-{figTag}.eps", and initializes it.
    Returns the {PSStream}.
    
    The nominal client plot window is {bbox}, plus a 
    safety margin {o->margin}. */

void finish_figure(PSStream *ps, PlotOptions *o);
  /* Finalizes the current figure written to file {ps}.
    Closes the file and returns NULL.
    Otherwise simply terminates the current page. */

void hatch_rectangle(PSStream *ps, interval_t B[], double hstep);
  /* Draws hatch-lines across the cell {B}, spaced {hstep} in client
    units, with the current pen. */

#endif
