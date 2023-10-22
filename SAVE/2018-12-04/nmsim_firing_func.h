#ifndef nmsim_firing_func_H
#define nmsim_firing_func_H

#define _GNU_SOURCE
#include <stdint.h>
 
/* Firing functions for the Galves-Löcherbach neuron model. */
/* Last edited on 2017-07-24 21:35:56 by stolfilocal */

typedef struct nmsim_firing_func_t nmsim_firing_func_t;
  /* A descriptor of of a firing function {Phi}. */
  
struct nmsim_firing_func_t *nmsim_firing_func_new
  ( char *class, 
    double V_M, 
    double D_M, 
    int32_t deg
  );
  /* Returns a descriptor for a firing function of the given {class},
    with midpoint potential {V_M} and slope
    {1/D_M} at that potential, and degree {deg}. The current classes are described in 
    {nmsim_firing_func_INFO}. */

void nmsim_firing_func_eval
  ( struct nmsim_firing_func_t *Phi,
    double V, 
    double *pr,
    double *dpr
  );
  /* Evaluates the firing function {Phi} at the argument {V} (which
    should have been multiplied by the firing modulator).
    
    If {dpr} is {NULL}, the derivative will not be computed. */

double nmsim_firing_func_min_V(struct nmsim_firing_func_t *Phi);
  /* A potential at which the probability of firing is negligible. */

double nmsim_firing_func_max_V(struct nmsim_firing_func_t *Phi);
  /* A potential at which the probability of firing is 1 for all purposes. */

nmsim_firing_func_t *nmsim_firing_func_gauss_new(double V_M, double D_M);
  /* Returns a descriptor for a Gaussian-integral firing
    function {Phi} with midpoint potential {V_M} and slope
    {1/D_M} at that potential. */

#define nmsim_firing_func_INFO \
  "??? The available firing functions are \"gauss\" (integral of Gaussian bell) ... ???"

#endif
