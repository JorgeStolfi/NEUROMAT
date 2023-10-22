#ifndef nmsim_neuron_net_H
#define nmsim_neuron_net_H
 
/* Neuron-level modeling of networks of Galves-Löcherbach neurons. */
/* Last edited on 2017-07-24 23:48:43 by stolfilocal */

#define _GNU_SOURCE
#include <stdint.h>

#include <nmsim_neuron_parms.h>
#include <nmsim_neuron_state.h>

/* !!! remove state from network description !!! */
  
typedef uint64_t nmsim_neuron_count_t;
  /* A count of neurons. */
  
typedef uint64_t nmsim_neuron_ix_t;
  /* Index of a neuron in a networks. */
  
typedef uint64_t nmsim_synapse_count_t;
  /* A count of synapses. */
  
typedef uint64_t nmsim_synapse_ix_t;
  /* Index of a synapse among all neurons or in a specific neuron. */

#define nmsim_neuron_net_MAX_NEURONS (((uint64_t)1) << 36)
  /* Max number of neurons in a network. */

#define nmsim_neuron_net_MAX_SYNAPSES (((uint64_t)1) << 36)
  /* Max number of synapses in a network. */

typedef struct nmsim_neuron_net_t
  { nmsim_neuron_count_t N;  /* Number of neurons. */
    /* These arrays are indexed with neuron index {0..N-1}: */
    nmsim_neuron_state_t *state; /* Current state of each neuron. */
    nmsim_neuron_parms_t **parms; /* Parameters of each neuron. */
    nmsim_neuron_count_t *deg_in;  /* Number of input synapses of each neuron. */
    nmsim_neuron_count_t *deg_ot;  /* Number of output synapses of each neuron. */
    /* These arrays are indexed with {i} in {0..N-1} and synapse index in {0..deg_in[i]-1}: */
    nmsim_neuron_ix_t **src_in; /* Indices of source neurons of input synapses. */
    double **wt_in;    /* Rest weights of input synapses. */
    /* These parameters are indexed with {i} in {0..N-1} and synapse index in {0..deg_ot[i]-1}: */
    nmsim_neuron_ix_t **dst_ot; /* Indices of destination neurons of output synapses. */
  } nmsim_neuron_net_t;
  /* Description of a network with {N} Galves-Loecherbach neurons,
    with specific synapses and neuron and synapse parameters.
    
    The pointers {parms[i]} and {parms[j]} for distinct neurons {i,j}
    may point to the same parameter record. */
    
/* NETWORK CREATION */

nmsim_neuron_net_t *nmsim_neuron_net_new(nmsim_neuron_count_t N);
  /* Allocates an {nmsim_neuron_net_t} structure with {N}
    neurons.  Initially each neuron has zero synapses,
    all parameter and synapse table pointers are {NULL},
    and all neuron states are undefined. */

void nmsim_neuron_net_create_synapses
  ( nmsim_neuron_net_t *net,
    nmsim_synapse_count_t M,
    double wt_avg,
    double wt_dev,
    double pr_neg
  );
  /* Creates {M} random synapses from neurons of {net}. Every pair of
    pre- and post-synaptic neurons is generated with the same
    probability, including connections from a neuron to itself and
    multiple connections between the same pair of neurons.
    
    The synapse lists {net.syn_in} and {net.syn_ot} are left unsorted.
    The logarithms of the absolute weights will have a normal
    distribution with mean {A = log(wt_avg)} and standard deviation {D =
    log((wt_avg+wt_dev)/wt_avg)}.
    
    A fraction {pr_neg} of the neurons will have its output weights negated. */

/* NEURON-LEVEL SIMULATION */  

double nmsim_neuron_net_tot_input(nmsim_neuron_net_t *net, bool_t X[], nmsim_neuron_ix_t i);
  /* Computes the total voltage increment for neuron {i} of network {net} accumulated
    during a simulation step, given the firing indicators {X[0..net.N-1]} of all neurons 
    during that step. */

void nmsim_neuron_net_evolve
  ( void
  );
  /* Simulates the evolution of a network with one or more homogeneous
    populations of neurons with stochastic states {net.pso[k]} and
    neuron parameters {parm} from some discrete time {t} to time
    {t+1}.

    In addition to synaptic inputs, each neuron {i} in the 
    population {k} is also assumed to have received, between times
    {t} and {t+1}, an external input signal {I[i,t]},
    which is a normally distributed variable with mean
    and deviation given by {net.exin(k,t,&I_avg,I_dev)}

    ??? Should the input signal be modulated by the input modulator
    {G[i,t]} too? ???  */

#endif

