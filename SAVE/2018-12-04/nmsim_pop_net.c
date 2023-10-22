/* See {nmsim_pop_net.h} */
/* Last edited on 2017-07-24 21:36:14 by stolfilocal */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <affirm.h>
#include <bool.h>

#include <nmsim_neuron_parms.h>
#include <nmsim_cohort.h>
#include <nmsim_pop.h>

#include <nmsim_pop_net.h>

void nmsim_pop_net_mf_evolve
  (
  )
  {

    /* Compute the total input that each neuron {i} 
    in the population received from all its input synapses betwen time {t}
    and time {t+1} is a normal random variable
    with mean {DV_avg} and deviation {DV_dev}.  This potential
    increment includes the pulses {X[j,t]} of each input neuron {j}
    in the network, modulated by its output modulator {H[j,t]} and
    by the rested synaptic gain {wfix[j,i]}; but not
    yet modulated by the input gain {G[i,t]} of neuron {i}. */


    /*  to the state of the cohort of 
    neurons with age 0 (those that just fired).
    Namely, the fields {cs.V,cs.G,cs.H} are set
    to the reset values {V_R,G_R,H_R}, respectively.
    The fraction {cs.eta} is set to {rho}. */
  }

