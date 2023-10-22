#ifndef nmsim_neuron_state_H
#define nmsim_neuron_state_H
 
/* State of a Galves-Löcherbach neuron. */
/* Last edited on 2017-07-24 21:39:15 by stolfilocal */

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>

#include <nmsim_neuron_parms.h>

#define nmsim_neuron_state_FILE_TYPE "nmsim_neuron_state"
    
#define nmsim_neuron_state_VERSION "2017-07-21"

typedef struct nmsim_neuron_state_t
  { double V; /* Membrane potential. */
    uint64_t age; /* Number of steps since last firing. */
  } nmsim_neuron_state_t;

double nmsim_neuron_state_leakage(nmsim_neuron_parms_t *parms, uint64_t age);
double nmsim_neuron_state_input_gain(nmsim_neuron_parms_t *parms, uint64_t age);
double nmsim_neuron_state_output_gain(nmsim_neuron_parms_t *parms, uint64_t age);
  /* Return the leakage modulator {M}, input modulator {G}, and output modulator {H} of a 
    neuron with parameters {parms} and given {age}. */

void nmsim_neuron_state_throw(nmsim_neuron_state_t *state, nmsim_neuron_parms_t *parms, double rho);
  /* Sets {*state} to a random state compatible with {parms}.
    Assumes that the average activity of the neuron in the past was {rho},
    and therefore sets the firing age to about {1/rho}. */ 

#endif
