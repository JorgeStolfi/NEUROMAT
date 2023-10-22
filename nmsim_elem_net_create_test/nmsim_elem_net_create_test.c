#define PROG_NAME "nmsim_elem_net_create_test"
#define PROG_DESC "creates various test networks for {nmsim_elem_net_simulate}"
#define PROG_VERS "1.0"

/* Last edited on 2021-01-09 16:53:42 by jstolfi */ 

#define PROG_COPYRIGHT \
  "Copyright © 2020  State University of Campinas (UNICAMP)"
  
#define PROG_AUTH \
  "Created by J. Stolfi on 2020-12-24"
  
#define PROG_HIST
  
#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program writes to standard output" \
  " a simple element-level network, suitable for {nmsim_elem_net_simulate}. \n" \
  "\n" \
  "  " nmsim_elem_net_band_make_INFO "\n" \
  "\n" \
  "OPTIONS\n" \
  "  -numLayers {numLayers} \n" \
  "    This mandatory command line argument specifies the" \
  " the number of layers in each independent sub-network.\n" \
  "\n" \
  "  -bandWidth {bandWidth} \n" \
  "    This mandatory command line argument specifies the" \
  " number of neurons in each layer of each independent sub-network.\n" \
  "\n" \
  "  -numBands {numBands} \n" \
  "    This mandatory command line argument specifies the" \
  " number of independent copies (sub-networks) in the output .\n" \
  "\n" \
  "  -closed {closed} \n" \
  "    This mandatory Boolean command line argument specifies" \
  " whether there should be " \
  " an extra layer of synapses from the last layer back to" \
  " the first one. The value can be \"0\" or \"F\" for" \
  " false, \"1\", or \"T\" for true.\n" \
  "\n" \
  "  -phiType {phiType} \n" \
  "    This mandatory command line argument specifies the" \
  " type of the firing function {Phi}, encoded as a single letter:\n" \
  "\n" \
  "" nmsim_firing_func_class_INFO "\n" \
  "\n" \
  "  -V_tau {V_tau} \n" \
  "    This mandatory command line argument specifies the" \
  " characteristic time for recharging to the resting potential.\n" \
  "\n" \
  "  -M_R {M_R} \n" \
  "    This optional command line argument specifies the value of the" \
  " recharge modulator {M} just after firing.  If not specified, 1 is" \
  " assumed (no modulation).\n" \
  "\n" \
  "  -M_tau {M_tau} \n" \
  "    This optional command line argument specifies the characteristic" \
  " time for the restoration of the recharge modulator {M} towards the" \
  " resting value of 1.  It is only meaningful if {M_R} is not 1.  If not" \
  " specified, 0 is assumed (immediate restoration in the next step).\n" \
  "\n" \
  "  -WAvg {WMin} {WMax} \n" \
  "    This command line argument specifies the range of values for the" \
  " average {W_avg}  of the strength of the synapses into each neuron. It is" \
  " mandatory if the network has any synapses.  Each sub-network will have a" \
  " value of {W_avg} that interpoates between {WMin/N} and {WMax/N} where {N} is" \
  " the number of inputs of the neurons.  If there is only one" \
  " sub-network, the mean of the two values is used.\n" \
  "\n" \
  "  -WDev {WDev} \n" \
  "    This optionalcommand line argument specifies the relative standard deviation for" \
  " the synaptic strengths.  The strengths of the synapses into each neuron will" \
  " be drawn from a log-normal distribution with a mean {W_avg} that is specific" \
  " to each sub-network, and standard deviation {W_dev = W_avg*WDev}.  If not" \
  " specified, zero is assumed, meaning that all synapses will have the same weight {W_avg}.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  nmsim_elem_net_simulate(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2020-12-24 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2020-12-24 J. Stolfi: created.\n" \
  "  2021-01-07 J. Stolfi: changed to general multiple bands structure.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " PROG_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define stringify(x) strngf(x)
#define strngf(x) #x

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <argparser.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <affirm.h>

#include <nmsim_basic.h>
#include <nmsim_class_net.h>
#include <nmsim_group_net.h>
#include <nmsim_elem_net.h>
#include <nmsim_elem_net_band.h>

typedef struct nenct_options_t
  { int32_t numBands;  /* Number of independent sub-networks (bands). */
    int32_t numLayers;  /* Number of layes (neuron groups). */
    int32_t bandWidth;  /* Number of neurons in each layr of each sub-network. */
    bool_t closed;      /* If true, the last layer feeds back on the first one. */
    /* Neuron paramters: */
    char phiType;       /* Type of the firing function {Phi}. */
    double V_tau;       /* Characteristic time for recharging to the resting potential. */
    double M_R;         /* Value of the recharge modulator {M} just after firing. */
    double M_tau;       /* Characteristic time for the restoration of{M} towards 1. */
    /* Synapse parameters: */
    double WMin;        /* Min value for the average total synaptic input strength per neuron. */
    double WMax;        /* Max value for the average total synaptic input strength per neuron. */
    double WDev;        /* Relative deviation of synaptic input strengths. */
  } nenct_options_t;
  /* Arguments from command line. */
 
nmsim_elem_net_t *nenct_create_net(nenct_options_t *o);
  /* Builds and returns the test network specified by {o}. */

nmsim_class_neuron_t *nenct_make_class_neuron(char phiType, double V_tau, double M_R, double M_tau);
  /* Creates a neuron class with the given parameters and default values for the
    others. */
    
nmsim_class_synapse_t *nenct_make_class_synapse(double WMin, double WMax);
  /* Creates a synapse class with the given parameters and default values for the
    others.  The {W_avg} is the average of {WMin} and {WMax}, and the 
    deviation is that of a uniform variate over that range. */
    
nenct_options_t *nenct_parse_options(int argc, char **argv);
  /* Parses the command line arguments. */

int main(int argc, char **argv);

/* IMPLEMENTATIONS: */

#define timeStep_DEFAULT 1.0
  /* Timestep to assume when converting between {_mu} parameters and {_tau}
    parameters. */

int main(int argc, char **argv)
  { 
    nenct_options_t *o = nenct_parse_options(argc, argv);
    nmsim_elem_net_t *enet = nenct_create_net(o);
    nmsim_elem_net_write(stdout, enet, timeStep_DEFAULT);
    return 0;
  }
  
nmsim_elem_net_t *nenct_create_net(nenct_options_t *o)
  {
    if (o->closed || (o->numLayers > 1))
      { /* Network has synapses: */
        demand(! isnan(o->WMin), "must define the min synapse strength with \"-WAvg\"");
        demand(! isnan(o->WMax), "must define the max synapse strength with \"-WAvg\"");
      }

    /* Create the class-level network: */
    nmsim_class_neuron_count_t  nnc = 1;  /* Number of neuron classes. */
    nmsim_class_synapse_count_t nsc = 1;  /* Number of synapse classes. */

    nmsim_class_net_t *cnet = nmsim_class_net_new(nnc, nsc);
    cnet->nclass[0] = nenct_make_class_neuron(o->phiType, o->V_tau, o->M_R, o->M_tau);
    cnet->sclass[0] = nenct_make_class_synapse(o->WMin, o->WMax);
    
    nmsim_elem_net_t *enet = 
      nmsim_elem_net_band_make
        ( cnet, o->numBands, o->numLayers, o->bandWidth, o->closed, o->WMin, o->WMax, o->WDev );
    return enet;
  }

nmsim_class_neuron_t *nenct_make_class_neuron(char phiType, double V_tau, double M_R, double M_tau)
  {
    double V_B = 0; /* As in Nilton's nets. */
    double V_R = 0; /* As in Nilton's nets. */
    double c_B = nmsim_basic_mu_from_tau(V_tau, timeStep_DEFAULT);
    double M_mu = nmsim_basic_mu_from_tau(M_tau, timeStep_DEFAULT);
    double H_R = 1.0;
    double H_mu = 0.0; 
    double V_M = 20.0;
    double V_D =  5.0;
    nmsim_firing_func_t Phi = nmsim_firing_func_make(phiType, V_M, V_D);
    nmsim_class_neuron_t *nclass = nmsim_class_neuron_new(V_B, V_R, c_B, M_R, M_mu, H_R, H_mu, &Phi);
    return nclass;
  }
    
nmsim_class_synapse_t *nenct_make_class_synapse(double WMin, double WMax)
  { double W_avg = (WMin + WMax)/2;
    double W_dev = (WMax - WMin)/sqrt(12);
    nmsim_class_synapse_t *sclass = nmsim_class_synapse_new(W_avg, W_dev);
    return sclass;
  }
    
nenct_options_t *nenct_parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    /* Allocate the command line argument record: */
    nenct_options_t *o = notnull(malloc(sizeof(nenct_options_t)), "no mem");

    /* Parse keyword parameters: */

    argparser_get_keyword(pp, "-numBands");
    o->numBands = (int32_t)argparser_get_next_int(pp, 1, nmsim_elem_neuron_count_MAX);

    argparser_get_keyword(pp, "-numLayers");
    o->numLayers = (int32_t)argparser_get_next_int(pp, 1, nmsim_elem_neuron_count_MAX);

    argparser_get_keyword(pp, "-bandWidth");
    o->bandWidth = (int32_t)argparser_get_next_int(pp, 1, nmsim_elem_neuron_count_MAX);

    argparser_get_keyword(pp, "-closed");
    o->closed = argparser_get_next_bool(pp);
    
    argparser_get_keyword(pp, "-phiType");
    char *str = argparser_get_next_non_keyword(pp);
    if (strlen(str) != 1) { argparser_error(pp, "invalid {phiType}"); }
    o->phiType = str[0];
    
    argparser_get_keyword(pp, "-V_tau");
    o->V_tau = argparser_get_next_double(pp, 0.0, 10000.0);
     
    if (argparser_keyword_present(pp, "-M_R"))
      { o->M_R = argparser_get_next_double(pp, 0.0, 2.0);  }
    else
      { o->M_R = 1.0; }
     
    if (argparser_keyword_present(pp, "-M_tau"))
      { o->M_tau = argparser_get_next_double(pp, 0.0, 10000.0);  }
    else
      { o->M_tau = 0.0; }
    
    if (argparser_keyword_present(pp, "-WAvg"))
      { o->WMin = argparser_get_next_double(pp, -50.0, +50.0);
        o->WMax = argparser_get_next_double(pp, o->WMin, +50.0);
      }
    else
      { o->WMin = NAN; o->WMax = NAN; }
    
    if (argparser_keyword_present(pp, "-WDev"))
      { o->WDev = argparser_get_next_double(pp, 0.0, 5.0); }
    else
      { o->WDev = 0.0; }

    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }
