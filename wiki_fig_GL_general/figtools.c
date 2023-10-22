/* See figtools.h */
/* Last edited on 2019-10-19 03:17:24 by jstolfi */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

#include <affirm.h>
#include <jsstring.h>
#include <pswr.h>
#include <frgb.h>
#include <interval.h>

#include <figtools.h>

/* INTERNAL PROTOTYPES */

/* IMPLEMENTATIONS */

PSStream *start_figure
  ( PlotOptions *o,
    char *figTag,
    interval_t bbox[2],    /* Client plot window. */
    double scale           /* Plot scale (mm per client unit). */
  )
  { 
    fprintf(stderr, "=== %s =========================\n", figTag);
    
    double xsz = HI(bbox[0]) - LO(bbox[0]);
    double ysz = HI(bbox[1]) - LO(bbox[1]);

    /* Figure size in pt: */
    double hsz = scale * xsz * mm;
    double vsz = scale * ysz * mm;

    /* Figure margin width in pt: */
    double hmg = o->margin[0] * mm;
    double vmg = o->margin[1] * mm;
    
    char *psname = txtcat(o->outName, "-");
    PSStream *ps = pswr_new_stream(psname, NULL, TRUE, NULL, NULL, TRUE, 640.0, 480.0);
    affirm(ps != NULL, "no output stream");
    double htot = hsz + 2*hmg;
    double vtot = vsz + 2*vmg;
    pswr_set_canvas_size(ps, htot, vtot);
    pswr_new_canvas(ps, figTag);
    pswr_set_canvas_layout(ps, hsz, vsz, FALSE, hmg, vmg, 0, 1, 1);
    pswr_new_picture(ps, LO(bbox[0]), HI(bbox[0]), LO(bbox[1]), HI(bbox[1]));
    pswr_set_pen(ps, 0,0,0, 0.15, 0,0);
    return ps;
  }
  
void finish_figure(PSStream *ps, PlotOptions *o)
  { pswr_close_stream(ps); }

void hatch_rectangle(PSStream *ps, interval_t B[], double hstep)
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
        pswr_segment(ps, 
            xm - ta*nx, ym - ta*ny,
            xm + tb*nx, ym + tb*ny
          );
      }
  }
