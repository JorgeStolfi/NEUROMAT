#define PROG_NAME "nmsim_multipop"
#define PROG_DESC "topological slicing of triangle meshes for 3D printing"
#define PROG_VERS "1.0"

#define nmsim_multipop_C_COPYRIGHT \
  "Copyright © 2016 by University of São Paulo (USP) and State University of Campinas (UNICAMP)"

/* Last edited on 2016-09-28 17:14:25 by stolfilocal */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -netFile {NET_FILE} \\\n" \
  "    -timeStep {TIME_STEP} \\\n" \
  "    -nSteps {N_STEPS} \\\n" \
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
  "  The program reads a statistical description of a large network of Galves-Loecherbach (GN) neurons, and simulates its behavior for a specified number of time steps.\n" \
  "\n" \
  "  "  nmsim_net_read_INFO "\n" \
  "\n" \
  "  ???.\n" \
  "\n" \
  "OUTPUT FILES\n" \
  "  " nmsim_multipop_log_INFO "\n" \
  "\n" \
  "  " nmsim_multipop_phase_data_INFO "\n" \
  "\n" \
  "OPTIONS\n" \
  "  -netFile {NET_FILE}\n" \
  "    This mandatory argument specifies the name of the input file that describes the network.\n" \
  "\n" \
  "  -timeStep {TIME_STEP} \n" \
  "    This mandatory argument is the time step (in seconds) to be assumed in the simulation.\n" \
  "\n" \
  "  -nSteps {N_STEPS}\n" \
  "    This mandatory argument specifies the number of time steps to execute in the simulation.\n" \
  "\n" \
  "  -outPrefix {OUT_PREFIX}\n" \
  "    This mandatory argument specifies the common prefix for all output files.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  bolognac(1), prosciuttox(1), pancettaz(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2016-08-11 by Jorge Stolfi.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2016-08-11 Created by J.Stolfi.\n"             \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " nmsim_multipop_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define stringify(x) strngf(x)
#define strngf(x) #x

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <vec.h>
#include <argparser.h>
#include <affirm.h>

#include <libnmsim.h>

#include <nmsim_multipop_planes.h>
#include <nmsim_multipop_stats.h>
#include <nmsim_multipop_slicer.h>
#include <nmsim_multipop_closer.h>
/* #include <nmsim_multipop_xml.h> */
/* #include <nmsim_multipop_graphics.h> */

typedef struct nmsim_multipop_options_t
  {
    /* Input: */
    char netFile;         /* Name of file with network description. */
    
    /* Simulation parameters: */
    double timeStep;      /* Time step (seconds). */
    int64_t nSteps;       /* Number of tim steps to simulate. */

    /* Output files: */
    char *outPrefix;  /* Prefix for output file names. */

  } nmsim_multipop_options_t;

/* INTERNAL PROTOTYPES */

int main(int argc, char **argv);
  /* Main program. */

nmsim_multipop_options_t *nmsim_multipop_parse_options(int argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

/* IMPLEMENTATIONS */

int main(int argc, char **argv)
  {
    /* Parse the command line options: */
    nmsim_multipop_options_t *o = nmsim_multipop_parse_options(argc, argv);

    /* Get the input triangle mesh {mesh}, quantizing to even integers: */
     mesh = stmesh_read_STL(o->modelFile, o->binary, o->eps, o->nfGuess, TRUE, o->preSorted);

    /* Get the Z-coordinates of the slicing planes, in increasing order: */
    i3_t minP, maxP; /* Bounding box of mesh. */
    stmesh_get_bounding_box(mesh, &minP, &maxP);
    bool_t uniform = ((! isnan(o->deltaZ)) && (o->deltaZ != 0));
    int_vec_t planeZ; /* The {Z}-coordinates of the slicing planes, quantized by {eps}. */
    if (uniform)
      { planeZ = nmsim_multipop_planes_get_uniform(o->startZ, o->deltaZ, o->eps, minP.c[2], maxP.c[2]); }
    else
      { planeZ = nmsim_multipop_planes_get_adaptive(o->planeZFile, o->eps, minP.c[2], maxP.c[2]); }

    /* Process it: */
    nmsim_multipop_stats_t st;
    nmsim_multipop_slicer_slice(mesh, o->preSorted, &planeZ, uniform, o->slicer, o->closer, o->outPrefix, &st);

    /* Show statistics: */
    nmsim_multipop_stats_print(stderr, &st);

    fprintf(stderr, "done.\n");
    return 0;
  }

/* ARGUMENT PARSING */
  
#define nmsim_multipop_timeStep_MIN (1.0e-6) 
  /* Minimum simulation time step (seconds). */
   
#define nmsim_multipop_timeStep_MAX (1.0e+2)
  /* Maximum simulation time step (seconds). */
  
#define nmsim_multipop_nSteps_MIN (0) 
  /* Minimum simulation step count. */
   
#define nmsim_multipop_nSteps_MAX (1024*1024*1024*1024-1)
  /* Maximum simulation step count. */
 
nmsim_multipop_options_t *nmsim_multipop_parse_options(int argc, char** argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    nmsim_multipop_options_t *o = notnull(malloc(sizeof(nmsim_multipop_options_t)), "no mem"); 
    
    /* Parse keyword parameters: */
    
  "    -outPrefix {OUT_PREFIX} \\\n" \


    /* Input network description file: */
    argparser_get_keyword(pp, "-netFile");
    o->netFile = argparser_get_next_non_keyword(pp);
    
    /* Time step: */
    argparser_get_keyword(pp, "-timeStep");
    o->timeStep = argparser_get_next_double(pp, nmsim_multipop_timeStep_MIN, nmsim_multipop_timeStep_MAX);
    
    /* Number of steps: */
    argparser_get_keyword(pp, "-nSteps");
    o->nSteps = argparser_get_next_int(pp, nmsim_multipop_nSteps_MIN, nmsim_multipop_nSteps_MAX);

    /* Outpu file prefix: */
    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next_non_keyword(pp);
    
    /* Parse positional arguments: */
    argparser_skip_parsed(pp);

    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }

