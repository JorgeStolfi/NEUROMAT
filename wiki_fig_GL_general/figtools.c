/* See figtools.h */
/* Last edited on 2023-12-08 08:39:19 by stolfi */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

#include <affirm.h>
#include <jsstring.h>
#include <epswr.h>
#include <frgb.h>
#include <interval.h>

#include <figtools.h>

/* INTERNAL PROTOTYPES */

/* IMPLEMENTATIONS */

epswr_figure_t *figtools_start_figure
  ( figtools_options_t *o,
    char *figTag,
    interval_t bbox[2]     /* Client plot window (client units). */
  )
  { 
    fprintf(stderr, "=== %s =========================\n", figTag);
    
    double xsz = HI(bbox[0]) - LO(bbox[0]); /* Plot window width (client units) */
    double ysz = HI(bbox[1]) - LO(bbox[1]); /* Plot window height (client units) */

    /* Figure size in pt: */
    double hsz = o->scale * xsz * epswr_pt_per_mm;
    double vsz = o->scale * ysz * epswr_pt_per_mm;

    /* Figure margin width in pt: */
    double hmg = o->margin[0] * epswr_pt_per_mm;
    double vmg = o->margin[1] * epswr_pt_per_mm;

    char *fname = txtcat(o->outName, ".eps");
    epswr_figure_t *eps = epswr_new_named_figure(NULL, NULL, o->outName, -1, figTag, hsz, vsz, hmg, hmg, vmg, vmg, TRUE); 
    epswr_set_client_window(eps, LO(bbox[0]), HI(bbox[0]), LO(bbox[1]), HI(bbox[1]));
    
    epswr_set_pen(eps, 0,0,0, 0.15, 0,0);
    return eps;
  }
  
void figtools_finish_figure(epswr_figure_t *eps, figtools_options_t *o)
  { epswr_end_figure(eps); }

void figtools_hatch_rectangle(epswr_figure_t *eps, interval_t B[], double hstep)
  {
    int i;
    double dx = HI(B[0]) - LO(B[0]);
    double dy = HI(B[1]) - LO(B[1]);
    double nx = -dy, ny = dx;
    int nh = (int)ceil(hypot(dx,dy)/hstep);
    fprintf(stderr, 
      "hatch hstep = %.5f d = (%.5f, %.5f) nh = %d\n",
      hstep, dx, dy, nh
    );
    if (nh < 3) { nh = 3; }
    for (i = 1; i < nh; i++)
      { double t = ((double) i)/((double)nh);
        double xm = LO(B[0]) + t*dx;
        double ym = LO(B[1]) + t*dy;
        double lim;
        /* Find {ta} st. {(xm,ym) - ta*(nx,ny)} is on low-right border */
        double ta = 4.0;
        lim = xm - HI(B[0]); if (ta*nx < lim) { ta = lim/nx; }
        lim = ym - LO(B[1]); if (ta*ny > lim) { ta = lim/ny; }
        /* Find {tb} st. {(xm,ym) + ta*(nx,ny)} is on up-left border */
        double tb = 4.0;
        lim = LO(B[0]) - xm; if (tb*nx < lim) { tb = lim/nx; }
        lim = HI(B[1]) - ym; if (tb*ny > lim) { tb = lim/ny; }
        /* Plot segment: */
        epswr_segment(eps, 
            xm - ta*nx, ym - ta*ny,
            xm + tb*nx, ym + tb*ny
          );
      }
  }
