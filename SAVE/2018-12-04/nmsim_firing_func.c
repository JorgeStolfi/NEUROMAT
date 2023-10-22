/* See {nmsim_firing_func.h} */
/* Last edited on 2017-07-24 21:51:52 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#include <affirm.h>
#include <bool.h>

#include <nmsim_firing_func.h>

typedef void nmsim_firing_func_eval_proc_t
  ( double V, 
    double V_M, 
    double D_M,
    double deg,
    double *pr,
    double *dpr
  );
  /* Type of a procedure that implements a certain class
    of firing functions. 
  
    The call {eval(V,V_M,D_M,deg,&pr,&dpr)} should return the firing
    probability {*pr=Phi(V)} and its derivative {*dpr = \partial_V
    Phi(V)} for a given membrane potential {V} (which should have been
    mutiplied by the firing modulator).
    
    The procedure should return {*pr = 1/2} and {*dpr = D_M}, for any
    {deg}, when {V = V_M}.
    
    If {dpr} is {NULL}, the derivative will not be computed.
    
    The {deg} parameter depends on the class, and may be ignored. */

nmsim_firing_func_t *nmsim_firing_func_new_gen
  ( char *desc, 
    nmsim_firing_func_eval_proc_t *eval, 
    double deg, 
    double V_M, 
    double D_M
  );
  /* Creates a new {nmsim_firing_func_t} record with given fields. */ 

typedef struct nmsim_firing_func_t
  { nmsim_firing_func_eval_proc_t *eval;            /* Evaluates {Phi(V)}. */
    double V_M;      /* Midpoint reference potential. */
    double D_M;      /* The slope at the midpoint. */
    double deg;      /* A shape parameter. */
    char *desc;      /* Description of function. */
  } nmsim_firing_func_t;
  /* Description of a firing function {Phi}, namely the member of 
    the class defined by {eval} selected by the parameters
    {V_M,D_M,deg}.
    
    The {desc} field is a string that contains a human-readable
    description of the firing function, including the values of
    the parameters {V_M,D_M,deg} when relevant. */  
    
/* GENERIC PROCEDURES */

struct nmsim_firing_func_t *nmsim_firing_func_new
  ( char *class, 
    double V_M,
    double D_M, 
    int32_t deg
  )
  {
    nmsim_firing_func_t *Phi = NULL;
    if (strcmp(class, "Gauss") == 0)
      { demand(deg == 0, "invalid degree for Gauss Phi");
        Phi = nmsim_firing_func_gauss_new(V_M, D_M);
      }
    else
      { demand(FALSE, "invalid class"); }
    return Phi;
  }

void nmsim_firing_func_eval
  ( nmsim_firing_func_t *Phi,
    double V, 
    double *pr,
    double *dpr
  )
  {
    Phi->eval(V, Phi->V_M, Phi->D_M, Phi->deg, pr, dpr);
  }
 
/* SPECIFIC FIRING FUNCTIONS */    

void nmsim_firing_func_gauss_eval
  ( double V, 
    double V_M,
    double D_M,
    double deg,
    double *pr,
    double *dpr
  );

nmsim_firing_func_t *nmsim_firing_func_gauss_new(double V_M, double D_M)
  {
    char *desc = NULL;
    asprintf(&desc, "PhiGauss[V_M=%.3f,D_M=%.3f]", V_M, D_M);
    nmsim_firing_func_eval_proc_t *eval = &nmsim_firing_func_gauss_eval;
    nmsim_firing_func_t *Phi = nmsim_firing_func_new_gen(desc, eval, 0, V_M, D_M);
    return Phi;
  }

#define nmsim_firing_func_gauss_Z_SAT (6.0)
 /* Value of {z} parameter for which the Gaussian Phi saturates. */

void nmsim_firing_func_gauss_eval
  ( double V, 
    double V_M,
    double D_M,
    double deg,
    double *pr,
    double *dpr
  )
  {
    assert(deg == 0);
    double z = (V - V_M)/D_M;
    if (z < -nmsim_firing_func_gauss_Z_SAT)
      { (*pr) = 0.0; (*dpr) = 0.0; }
    else if (z > +nmsim_firing_func_gauss_Z_SAT)
      { (*pr) = 1.0; (*dpr) = 0.0; }
    else
      { /* !!! FIX THIS: !!! */
        (*pr) = 0.5; (*dpr) = 0.0;
      }
  }

nmsim_firing_func_t *nmsim_firing_func_new_gen
  ( char *desc, 
    nmsim_firing_func_eval_proc_t *eval, 
    double deg, 
    double V_M, 
    double D_M
  )
  {
    nmsim_firing_func_t *Phi = notnull(malloc(sizeof(nmsim_firing_func_t)), "no mem");
    (*Phi) = (nmsim_firing_func_t){ .eval = eval, .V_M = V_M, .D_M = D_M, .deg = deg, .desc = desc };
    return Phi;
  }
