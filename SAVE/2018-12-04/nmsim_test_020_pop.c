#define PROG_NAME "nmsim_test_010_neuron"
#define PROG_DESC "basic tests of {limnmism} neuron-level procedures"
#define PROG_VERS "1.0"

/* Last edited on 2018-03-04 23:00:01 by stolfilocal */ 

#define PROG_COPYRIGHT \
  "Copyright © 2017  State University of Campinas (UNICAMP)"
  
#define PROG_AUTH \
  "Created by J. Stolfi on 2017-07-21"
  
#define PROG_HIST
  
#define PROG_HELP \
  "  " PROG_NAME ""

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <argparser.h>
#include <jsfile.h>
#include <affirm.h>

#include <nmsim_firing_func.h>
#include <nmsim_neuron_parms.h>
#include <nmsim_neuron_state.h>

void nmsim_test_discrete_net(uint64_t N, uint64_t T);
  /* Tests a neuron net with {N} neurons for {T} steps. */

int main(int argc, char **argv);

/* IMPLEMENTATIONS: */

int main(int argc, char **argv)
  { 
    nmsim_test_discrete_net(1,1000);
    nmsim_test_discrete_net(100,10000);
    nmsim_test_discrete_net(10000,100000);
    nmsim_test_discrete_net(1000000,1000000);
    return 0;
  }
    
void nmsim_test_discrete_net(uint64_t N, uint64_t T)
  {
    fprintf(stderr, "--- testing discrete network ---\n");
    fprintf(stderr, "network size {N} = %lu neurons\n", N);
    fprintf(stderr, "simulation length {T} = %lu steps\n", T);
    bool_t debug = FALSE;
    
    /* Creates a random network: */
    nmsim_neuron_net_t *net = nmsim_neuron_net_new(N);
    nmsim_test_net_init(net);
    nmsim_test_net_simulate(net, T);
  }

void nmsim_test_net_init
  ( nmsim_neuron_net_t *net,
    nmsim_neuron_params_t *parms
  )
  {
    for (uint64_t i = 0; i < net->N; i++) 
      { nmsim_neuron_state_t *state = &(net->state[i]); 
        /* Set the parameter pointer: */
        net->parms[i] = parms;
        /* Creates synapses: */
        nmsim_test_net_create_in_synapses(net, deg_in_min, deg_in_max);
       
        double rho = 0.005; /* For now, assume firing every 100 steps. */
        nmsim_neuron_state_throw(state, parms, rho);
      }
    nmsim_test_net_collect_ot_synapses(deg_ot_min, deg_in_max, &(deg_in[i]), &(syn_in[i]), N);


  }
    
void nmsim_neuron_net_create_synapses
  ( nmsim_neuron_net_t *net,
    uint64_t deg_in_min, 
    uint64_t deg_in_max
  )
  {
    &(deg_in[i]), &(syn_in[i]), N
    unit64_t deg_in = uint64_abrandom(deg_in_min, deg_in_max);
    uint64_t *syn_in = notnull(malloc(deg_in*sizeof(uint64_t)_, "no mem");
    for (uint64_t k = 0; k < deg_in; k++) 
      { syn_in[k] = uint64_abrandom(0, N-1); }
      
    (*deg_inP) = deg_in;
    (*syn_inP) = syn_in;
  }
