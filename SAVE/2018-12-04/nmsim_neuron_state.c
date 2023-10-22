/* See {nmsim_neuron_state.h} */
/* Last edited on 2017-07-25 00:34:02 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <nget.h>
#include <fget.h>
#include <filefmt.h>
#include <jsmath.h>
#include <jsrandom.h>

#include <nmsim_firing_func.h>
#include <nmsim_neuron_state.h>
 
double nmsim_neuron_state_exp_dyn(double R, double mu, uint64_t age);
  /* The value of a modulator with reset-decay parameters {R,mu}
    for a neuron with specified firing {age}. */

void nmsim_neuron_state_throw(nmsim_neuron_state_t *state, nmsim_neuron_parms_t *parms, double rho) 
  { 
    demand((rho > 1.0e-8) && (rho <= 1.0), "invalid mean activity");
    /* Random potential: */
    double V_sat = nmsim_firing_func_max_V(parms->Phi);
    double V_min = fmin(parms->V_B, parms->V_R);
    double V_max = fmax(fmax(parms->V_B, parms->V_R), V_sat);
    state->V = dabrandom(V_min, V_max);

    /* Random age: */
    uint64_t age = (uint64_t)floor(1.0/rho + 0.5);
    state-> age = age;
  }

double nmsim_neuron_state_leakage(nmsim_neuron_parms_t *parms, uint64_t age)
  { return nmsim_neuron_state_exp_dyn(parms->M_R, parms->M_mu, age); }
  
double nmsim_neuron_state_input_gain(nmsim_neuron_parms_t *parms, uint64_t age)
  { return nmsim_neuron_state_exp_dyn(parms->G_R, parms->G_mu, age); }

double nmsim_neuron_state_output_gain(nmsim_neuron_parms_t *parms, uint64_t age)
  { return nmsim_neuron_state_exp_dyn(parms->H_R, parms->H_mu, age); }
  
double nmsim_neuron_state_exp_dyn(double R, double mu, uint64_t age)
  {
    /* Should use a table and have a better cutoff. */
    if ((mu = 1.0) || (R == 1.0) || (age > 1000000)) { return 1.0; }
    double P = 1.0 - (1.0 - R)*pow(mu, (double)age);
    return P;
  }
