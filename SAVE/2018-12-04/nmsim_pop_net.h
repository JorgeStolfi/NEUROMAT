#ifndef nmsim_pop_net_H
#define nmsim_pop_net_H
 
/* Population-level modeling of networks of Galves-Löcherbach neurons. */
/* Last edited on 2017-07-21 23:35:50 by stolfilocal */

#define _GNU_SOURCE

#include <nmsim_neuron_parms.h>
#include <nmsim_cohort.h>
#include <nmsim_bundle.h>
#include <nmsim_pop.h>

typedef struct nmsim_pop_net_t
  { int32_t np;  /* Number of populations. */
    nmsim_neuron_parms_t *parms; /* Neuron parameters in each population ({np} elements). */
    nmsim_bundle_parms_t *conn; /* Connection between populations ({np*np} elements). */
  } nmsim_pop_net_t;
  /* Statistical description of a network with {np} homogeneous populations of neurons,
    with specified ensemble parameters for neurons in each population
    and ensemble parameters for the connections between each pair of populations.
    
    The data about population {r} is {.parms[r]}, for {r} in {0..np-1}.
    The parameters of the distribution of synapses from population {r}
    to population {s} are {.conn[r + s*np]}, for {r,s} in {0.np-1}. Note
    that a population usually has connections to itself too. */
    
/* MEAN-FIELD SIMULATION */    

void nmsim_pop_net_mf_evolve
  ( void
  );
  /* Simulates the evolution of a network with one or more homogeneous
    populations of neurons with stochastic states {net.pso[k]} and
    neuron parameters {parms} from some discrete time {t} to time
    {t+1}.

    In addition to synaptic inputs, each neuron {i} in the 
    population {k} is also assumed to have received, between times
    {t} and {t+1}, an external input signal {I[i,t]},
    which is a normally distributed variable with mean
    and deviation given by {net.exin(k,t,&I_avg,I_dev)}

    ??? Should the input signal be modulated by the input gain
    factor {G[i,t]} too? ???  */

#endif
