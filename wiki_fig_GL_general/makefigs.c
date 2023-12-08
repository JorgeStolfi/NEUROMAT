#define PROG_NAME "makefigs"
#define PROG_DESC "generate figures for the Wikipedia GL neuron model - general version"
#define PROG_VERS "1.0"

/* Copyright © 2004 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */
/* Last edited on 2023-12-08 10:38:53 by stolfi */

/* Created may/2004 by Jorge Stolfi, UNICAMP */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <values.h>
#include <assert.h>

#include <epswr.h>
#include <figtools.h>

/* INTERNAL PROTOTYPES */

int main(int argc, char **argv);

void bitstreams(figtools_options_t *o);
  /* Write the streams of bits figure. */

figtools_options_t *parse_options(int argc, char **argv);
void get_arg_double(double *varp, int *argnp, int argc, char **argv, char *usage);
void get_arg_string(char **varp, int *argnp, int argc, char **argv, char *usage);
void arg_error(char *msg, char *arg, char *pname, char *usage);

/* IMPLEMENTATIONS */

int main(int argc, char **argv)
  {   
    fprintf(stderr, "Parsing options ...\n");
    figtools_options_t *o = parse_options(argc, argv);
    
    fprintf(stderr, "Plotting figure {bitstreams} ...\n");
    bitstreams(o);
    
    fprintf(stderr, "Done.\n");
    return 0;
  }

void bitstreams(figtools_options_t *o)
  {
    
    int32_t rows =  7; /* Neurons. */
    int32_t cols = 18; /* Times. */

    double dx = 4.0; /* Client units. */
    double dy = 3.5; /* Client units. */

    double pw = 0.2; /* Pen width (client units). */

    auto void stateframe(double x, double y, frgb_t *color, double mrg);
      /* Draws the frame for the symbol box at {x,y), with color {color},
        expanded by {mrg} on both axes. */
    
    auto void predframe(int32_t iN, int32_t tN, frgb_t *color, double dxy);
      /* Draws the frame for neuron {iN} that last fired
        at time {tN}, with color {color}, displaced by {dxy}
        on both axes. */
    
    interval_t bbox[2];
    bbox[0] = (interval_t){{ -4.0*dx, +(cols + 1.0)*dx }};
    bbox[1] = (interval_t){{ -1.2*dy, +(rows + 0.3)*dy }};
    
    char *ftag = "bitstreams";
    fprintf(stderr, "  Creating %s-%s.eps ...\n", o->outName, ftag);
    epswr_figure_t *eps = figtools_start_figure(o, ftag, bbox);
    
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
        for (int32_t t = -3; t <= cols; t++)
          { double x = t*dx;
            
            /* Select what symbol to show, or none, and its font: */
            double fontSize = dy*3.1*o->scale;
            char *font = "Courier";
            char *dig = NULL;
            if (t < 0) 
              { dig = "."; }
            else if ((t >= 0) && (t < cols)) 
              { int32_t k = i*cols + t;
                dig = "0";
                if (X[k] != 0) { dig = "1"; font = "CourierBold"; }
              }
            epswr_set_label_font(eps, font, fontSize);
              
            /* Decide the symbol's color and whether to draw a box around it: */
            frgb_t color = colorN;
            double boxmrg = NAN; /* Margin for box, or {NAN} if no box. */
            if (t == cols)
              { if (i == iA) { color = colorA; }
                if (i == iB) { color = colorB; }
                boxmrg = 0.0;
              }
            else if ((i == iA) && (t == tA)) 
              { color = colorA;
                boxmrg = 0.0;
              }
            else if ((i == iB) && (t == tB))
              { color = colorB;
                boxmrg = 0.0;
              }
              
            if (! isnan(boxmrg)) { stateframe(x, y, &color, 0.0); }
            
            if (dig != NULL)
              { /* Draw the symbol: */
                epswr_set_fill_color(eps, color.c[0], color.c[1], color.c[2]);
                epswr_label(eps, dig, dig, x,y, 0.0, TRUE, 0.5,0.5, TRUE, FALSE);
              }
            x = x + dx;
          }
        y = y + dy;
      }
    
    predframe(iA, tA, &colorA, -1.0);
    predframe(iB, tB, &colorB, +1.0);
    
    // epswr_rectangle(eps, 
    //   LO(cell[0]), HI(cell[0]), 
    //   LO(cell[1]), HI(cell[1]), 
    //   TRUE, TRUE
    // );
    
    fprintf(stderr, "  Finishing the figure ...\n");
    figtools_finish_figure(eps, o);
    
    return;
    
    /* INTERNAL IMPLEMENTATIONS */
    
    void stateframe(double x, double y, frgb_t *color, double boxmrg)
      { /* Draw the box around the symbol: */
        /* Coordinates of sides of digit box: */
        double xblo = x - 0.30*dx - boxmrg;
        double xbhi = x + 0.30*dx + boxmrg;
        double yblo = y - 0.40*dy - boxmrg;
        double ybhi = y + 0.40*dy + boxmrg;
        epswr_set_pen(eps, color->c[0], color->c[1], color->c[2], pw*o->scale,  0.0,0.0);
        epswr_set_fill_color(eps, color->c[0], color->c[1], color->c[2]);
        epswr_rectangle(eps, xblo,xbhi, yblo,ybhi, FALSE, TRUE);
      }

    void predframe(int32_t iN, int32_t tN, frgb_t *color, double dxy)
      { 
        /* Coordinates of frame sides: */
        double xlo = (tN - 0.5)*dx;
        double xhi = (cols - 0.5)*dx + 0.3*dxy;
        
        double ylo = (0 - 0.7)*dy + 0.3*dxy;
        double yhi = (rows - 0.3)*dy + 0.3*dxy;
        
        double pw = 0.2; /* Pen width. */
        
        epswr_set_pen(eps, color->c[0],color->c[1],color->c[2], pw*o->scale,  0.0,0.0);
        epswr_set_fill_color(eps, color->c[0],color->c[1],color->c[2]);
        epswr_rectangle(eps, xlo, xhi, ylo, yhi, FALSE, TRUE);
        
        /* Coordinates of arrow: */
        double xa = xhi;
        double ya = (iN + 0.0)*dy;
        
        double xb = (cols + 0.0)*dx;
        double yb = ya;
        epswr_segment(eps, xa,ya, xb,yb);
        epswr_arrowhead(eps, xa,ya, xb,yb, 3.5, 5.5, 1.0, TRUE, TRUE);
      }
  }
  
#define ARG_ERROR(Msg,Arg) arg_error((Msg),(Arg),argv[0],usage)
  
#define GET_STRING(Var) get_arg_string(&(Var), &argn, argc, argv, usage)
#define GET_DOUBLE(Var) get_arg_double(&(Var), &argn, argc, argv, usage)

figtools_options_t *parse_options(int argc, char **argv)
  {
    figtools_options_t *o = (figtools_options_t *)malloc(sizeof(figtools_options_t));
    char* usage = "\n  [ -help ] [ -outName PREFIX ] [ -eps | -eps ]";
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
