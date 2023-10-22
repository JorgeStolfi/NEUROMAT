#define PROG_NAME "makefigs"
#define PROG_DESC "generate figures for the Wikipedia GL neuron model - general version"
#define PROG_VERS "1.0"

/* Copyright © 2004 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */
/* Last edited on 2019-10-19 06:23:28 by jstolfi */

/* Created may/2004 by Jorge Stolfi, UNICAMP */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <values.h>
#include <assert.h>

#include <pswr.h>
#include <figtools.h>

/* INTERNAL PROTOTYPES */

int main(int argc, char **argv);

void bitstreams(PlotOptions *o);
  /* Write the streams of bits figure. */

PlotOptions *parse_options(int argc, char **argv);
void get_arg_double(double *varp, int *argnp, int argc, char **argv, char *usage);
void get_arg_string(char **varp, int *argnp, int argc, char **argv, char *usage);
void arg_error(char *msg, char *arg, char *pname, char *usage);

/* IMPLEMENTATIONS */

int main(int argc, char **argv)
  {   
    fprintf(stderr, "Parsing options ...\n");
    PlotOptions *o = parse_options(argc, argv);
    
    fprintf(stderr, "Plotting figure {bitstreams} ...\n");
    bitstreams(o);
    
    fprintf(stderr, "Done.\n");
    return 0;
  }

void bitstreams(PlotOptions *o)
  {
    
    auto void predframe(int32_t iN, int32_t tN, frgb_t *color, double dxy);
      /* Draws the frame for neuron {iN} that last fired
        at time {tN}, with color {color}, displaced by {dxy}
        on both axes. */
    
    int32_t rows =  7; /* Neurons. */
    int32_t cols = 18; /* Times. */

    double dx = 4.0 * mm;
    double dy = 3.5 * mm;

    interval_t bbox[2];
    bbox[0] = (interval_t){{ -2.0*dx, +(cols + 1.0)*dx }};
    bbox[1] = (interval_t){{ -1.2*dy, +(rows + 0.3)*dy }};
    
    char *ftag = "bitstreams";
    fprintf(stderr, "  Creating %s-%s.eps ...\n", o->outName, ftag);
    PSStream *ps = start_figure(o, ftag, bbox, 1.0);
    
    /* Choose the neurons to show and their last firing times:*/
    int32_t iA = 1; int32_t tA = cols-6;
    int32_t iB = 4; int32_t tB = cols-4;
    
    /* Define the firing indicators and remember last firings of selected ones: */
    fprintf(stderr, "  Defining the firing history ...\n");
    int8_t X[rows*cols]; /* Firing indicators. */
    for (int32_t i = 0; i < rows; i++)
      { /* Write a row of 0's and 1's. */
        for (int32_t t = 0; t < cols; t++)
          { int32_t k = i*cols + t;
            double rn = 0.5*(1 + 0.6*sin(17.23*k) + 0.4*sin(33.22*t)); /* A random number. */
            X[k] = (int8_t)(rn > 0.6); /* A random bit. */
            if (i == iA)
              { if (t == tA) { X[k] = 1; }
                if (t > tA) { X[k] = 0;}
              }
            if (i == iB)
              { if (t == tB) { X[k] = 1; }
                if (t > tB) { X[k] = 0;}
              }
          }
      }
    
    /* Draw the firing history: */
    fprintf(stderr, "  Plotting the firing history ...\n");
    frgb_t colorN = (frgb_t){{ 0.000f, 0.000f, 0.000f }};
    frgb_t colorA = (frgb_t){{ 0.800f, 0.200f, 0.000f }};
    frgb_t colorB = (frgb_t){{ 0.000f, 0.200f, 1.000f }};
    for (int32_t i = 0; i < rows; i++)
      { /* Write a row of 0's and 1's. */
        double y = i*dy;
        for (int32_t t = -1; t <= cols; t++)
          { double x = t*dx;
            
            /* Select what symbol to show, or none, and its font: */
            char *font = "Courier";
            char *dig = NULL;
            if (t < 0) 
              { dig = "..."; }
            else if ((t >= 0) && (t < cols)) 
              { int32_t k = i*cols + t;
                dig = "0";
                if (X[k] != 0) { dig = "1"; font = "CourierBold"; }
              }
            pswr_set_label_font(ps, font, 18.0);
              
            /* Decide the symbol's color and whether to draw a box around it: */
            bool_t box = (t == cols); /* Should draw the digit box? */
            frgb_t color = colorN;
            double pw = 0.3; /* Pen width. */
            if ((i == iA) && ((t == tA) || (t == cols))) 
              { color = colorA; box = TRUE; pw = 0.7; }
            if ((i == iB) && ((t == tB) || (t == cols)))
              { color = colorB; box = TRUE; pw = 0.7; }
            
            pswr_set_pen(ps, color.c[0], color.c[1], color.c[2], pw,  0.0,0.0);
            pswr_set_fill_color(ps, color.c[0], color.c[1], color.c[2]);
            if (dig != NULL)
              { /* Draw the symbol: */
                pswr_fill_draw_label(ps, dig, x, y, 0.0, 0.5, 0.5, TRUE, FALSE);
              }
            if (box)
              { /* Draw the box around the symbol: */
                /* Coordinates of sides of digit box: */
                double xblo = x - 0.3*dx;
                double xbhi = x + 0.3*dx;
                double yblo = y - 0.4*dy;
                double ybhi = y + 0.4*dy;
                pswr_rectangle(ps, xblo,xbhi, yblo,ybhi, FALSE, TRUE);
              }
            x = x + dx;
          }
        y = y + dy;
      }
    
    predframe(iA, tA, &colorA, -1.0);
    predframe(iB, tB, &colorB, +1.0);
    
    // pswr_rectangle(ps, 
    //   LO(cell[0]), HI(cell[0]), 
    //   LO(cell[1]), HI(cell[1]), 
    //   TRUE, TRUE
    // );
    
    fprintf(stderr, "  Finishing the figure ...\n");
    finish_figure(ps, o);
    
    return;
    
    /* INTERNAL IMPLEMENTATIONS */
    
    void predframe(int32_t iN, int32_t tN, frgb_t *color, double dxy)
      { 
        /* Coordinates of frame sides: */
        double xlo = (tN - 0.3)*dx;
        double xhi = (cols - 0.7)*dx + dxy;
        
        double ylo = (0 - 0.5)*dy + 0.7*dxy;
        double yhi = (rows - 0.5)*dy + 0.7*dxy;
        
        double pw = 0.7; /* Pen width. */
        
        pswr_set_pen(ps, color->c[0],color->c[1],color->c[2], pw,  0.0,0.0);
        pswr_set_fill_color(ps, color->c[0],color->c[1],color->c[2]);
        pswr_rectangle(ps, xlo, xhi, ylo, yhi, FALSE, TRUE);
        
        /* Coordinates of arrow: */
        double xa = xhi;
        double ya = (iN + 0.0)*dy;
        
        double xb = (cols + 0.0)*dx;
        double yb = ya;
        pswr_segment(ps, xa,ya, xb,yb);
        pswr_arrowhead(ps, xa,ya, xb,yb, 2.0, 4.7, 1.0, TRUE, TRUE);
      }
  }
  
#define ARG_ERROR(Msg,Arg) arg_error((Msg),(Arg),argv[0],usage)
  
#define GET_STRING(Var) get_arg_string(&(Var), &argn, argc, argv, usage)
#define GET_DOUBLE(Var) get_arg_double(&(Var), &argn, argc, argv, usage)

PlotOptions *parse_options(int argc, char **argv)
  {
    PlotOptions *o = (PlotOptions *)malloc(sizeof(PlotOptions));
    char* usage = "\n  [ -help ] [ -outName PREFIX ] [ -eps | -ps ]";
    int argn;

    /* Defaults: */
    o->outName = "out/GL";
    o->color = TRUE;
    o->scale = 25.4;
    o->margin[0] = 1.0; 
    o->margin[1] = 1.0;

    argn = 1;

    /* Scan command line options. */
    while ((argn < argc) && (argv[argn][0] == '-') && (argv[argn][1] != '\0'))
      {
        char *key = argv[argn];
        if ((key[0] == '-') && (key[1] == '-') && (key[2] != '\0')) { key++; }
        if (strcmp(key, "-help") == 0)
          { fprintf(stderr, "usage: %s %s\n", argv[0], usage); exit(0); }
        else if (strcmp(key, "-outName") == 0)
          { GET_STRING(o->outName); }
        else if (strcmp(key, "-color") == 0)
          { o->color = TRUE; }
        else if ((strcmp(key, "-gray") == 0) || (strcmp(key, "-grey") == 0))
          { o->color = FALSE; }
        else if (strcmp(key, "-scale") == 0)
          { GET_DOUBLE(o->scale); }
        else if (strcmp(key, "-margin") == 0)
          { GET_DOUBLE(o->margin[0]);
            GET_DOUBLE(o->margin[1]);
          }
        else 
          { ARG_ERROR("unknown option", argv[argn]); }
        ++argn;
      }
    if (argn != argc) { ARG_ERROR("extraneous arguments", argv[argn]); }

    return o;
  }

void get_arg_string(char **varp, int *argnp, int argc, char **argv, char *usage)
   /*
     Stores the next command line argument (as a string) into "*varp" */
  {
    int argn = *argnp;
    if (argn+1 >= argc)
      { ARG_ERROR("missing arg value", argv[argn]); }
    (*varp) = argv[argn+1];
    (*argnp) = argn+1;
  }
  
void get_arg_double(double *varp, int *argnp, int argc, char **argv, char *usage)
   /*
     Stores the next command line argument (as a double) into "*varp" */
  {
    int argn = *argnp;
    char *end;
    if (argn+1 >= argc)
      { ARG_ERROR("missing arg value", argv[argn]); }
    (*varp) = strtod(argv[argn+1], &end);
    if ((*end) != '\0') 
      { ARG_ERROR("invalid numeric argument", argv[argn+1]); }
    (*argnp) = argn+1;
  }

void arg_error(char *msg, char *arg, char *pname, char *usage)
   /*
     Prints "msg", "arg", the "usage" message, and exits.
     Handy while parsing the command line arguments. */
  {
    fprintf(stderr, "%s %s\n", msg, arg);
    fprintf(stderr, "usage: %s %s\n", pname, usage);
    exit(1);
  }
