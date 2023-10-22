/* See {nmsim_cohort.h} */
/* Last edited on 2017-07-24 21:55:38 by stolfilocal */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>

#include <nmsim_firing_func.h>
#include <nmsim_neuron_parms.h>
#include <nmsim_cohort.h>

void nmsim_cohort_mf_state_set
  ( nmsim_cohort_mf_state_t *cs, 
    double V_avg,
    double V_dev,
    double M,
    double G,
    double H,
    double eta,
    nmsim_neuron_parms_t *parms
  )
  {
    double rho;
    nmsim_firing_func_eval(parms->Phi, G*V_avg, &rho, NULL);
    (*cs) = (nmsim_cohort_mf_state_t)
      { .V_avg = V_avg, .V_dev = V_dev,
        .M = M, .G = G, .H = H,
        .eta = eta, .rho = rho
      };
  }

void nmsim_cohort_mf_state_clear(nmsim_cohort_mf_state_t *cs)
  {
    (*cs) = (nmsim_cohort_mf_state_t)
      { .V_avg = 0.0, .V_dev = 0.0,
        .G = 1.0, .H = 1.0,
        .eta = 0.0, .rho = 0.0
      };
  }

void nmsim_cohort_mf_state_merge
  ( nmsim_cohort_mf_state_t *csa, 
    nmsim_cohort_mf_state_t *cst,
    nmsim_neuron_parms_t *parms
  )
  {
    double ea = csa->eta;
    double et = cst->eta;
    assert(ea >= 0);
    if (ea == 0) { return; }

    assert(et >= 0);
    if (et == 0)
      { /* Just replace: */
        (*cst) = (*csa); 
      }
    else
      { /* Get the main parameters: */
        double avga = csa->V_avg;
        double deva = csa->V_dev;
        double avgt = cst->V_avg;
        double devt = cst->V_dev;
        /* Total fraction of the two cohorts in population: */
        double eta_new = ea + et;
        /* Compute average potential of merged populations: */
        double avg_new = (ea*avga + et*avgt)/eta_new; 
        /* Compute the variance of the new population: */
        double da = avga - avg_new;
        double dt = avgt - avg_new;
        double dev_new = sqrt((ea*(da*da + deva*deva) + et*(dt*dt + devt*devt))/eta_new);
        /* Modulators (for now, to be deleted): */
        double M_new = (ea*(csa->M) + et*(cst->M))/eta_new;
        double G_new = (ea*(csa->G) + et*(cst->G))/eta_new;
        double H_new = (ea*(csa->H) + et*(cst->H))/eta_new;
        /* Store results: */
        nmsim_cohort_mf_state_set(cst, avg_new, dev_new, M_new, G_new, H_new, eta_new, parms);
      }        
  }

void nmsim_cohort_mf_evolve
  ( nmsim_cohort_mf_state_t *cso, 
    double DV_avg,
    double DV_dev,
    nmsim_neuron_parms_t *parms, 
    nmsim_cohort_mf_state_t *csn_fire,
    nmsim_cohort_mf_state_t *csn_fail
  )
  {
    /* Grab some variables for convenience: */
    double V_avg = cso->V_avg; /* Mean potential of cohort. */
    double V_dev = cso->V_dev; /* Potential deviation in cohort. */
    double M = cso->M;         /* Leakage modulator of cohort. */
    double G = cso->G;         /* Input gain of cohort. */
    double H = cso->H;         /* Output gain of cohort. */
    double eta = cso->eta;     /* Fraction of population in cohort. */
    double rho = cso->rho;     /* Fraction of cohort that is about to fire. */

    double c_B = parms->c_B;  /* Leakage at rest. */
    
    double M_mu = parms->M_mu;  /* Leakage modulator recovery factor. */
    double G_mu = parms->G_mu;  /* Input modulator recovery factor. */
    double H_mu = parms->H_mu;  /* Output modulator recovery factor. */

    /* Store the new state of the firing neurons in {csn_fire}: */
    double V_avg_fire = parms->V_R;
    double V_dev_fire = 0.0;
    nmsim_cohort_mf_state_set
      ( csn_fire,
        V_avg_fire, V_dev_fire, parms->M_R, parms->G_R, parms->H_R, eta*rho,
        parms
      );
   
    /* Compute the average potential of the non-firing neurons: */
    double c = c_B * M; /* Present leakage. */
    double V_avg_new = c*V_avg + DV_avg;
    /* Compute the deviation of the potential of those neurons: */
    /* ??? Can we ignore correlations here? ??? */
    double V_dev_new = hypot(c*V_dev, DV_dev);
    /* Compute the new leakage, input and output modulators of those neurons: */
    double M_new = 1.0 - M_mu*(1.0 - M);
    double G_new = 1.0 - G_mu*(1.0 - G);
    double H_new = 1.0 - H_mu*(1.0 - H);
    /* Store the new state of the non-firing neurons in {csn_fail}: */
    nmsim_cohort_mf_state_set
      ( csn_fail,
        V_avg_new, V_dev_new, M_new, G_new, H_new, eta*(1.0 - rho),
        parms
      );
  }

