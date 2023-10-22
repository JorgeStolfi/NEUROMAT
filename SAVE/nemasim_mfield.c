/* !!! Store {rho} in 

typedef struct nemasim_firing_function_t
  { nemasim_firing_proc_t *eval;  /* Procedure that evaluates {Phi(V)}. */
    char *name;                   /* Name of function. */
  } nemasim_firing_function_t;
  /* Parameters of a firing function. */    

typedef struct nemasim_neuron_parms_t
  { /* Parameters of potential dynamics: */
    double V_R;   /* Potential after firing. */
    double mu_V;  /* Potential decay factor. */
    /* Parameters of input gain dynamics: */
    double G_R;   /* Input gain after firing. */
    double mu_G;  /* Input gain recovery factor. */
    /* Parameters of output gain dynamics: */
    double H_R;   /* Output gain after firing. */
    double mu_H;  /* Output gain recovery factor. */
    /* Parameters of firing function: */
    nemasim_firing_function_t *Phi;
  } nemasim_neuron_parms_t;
  /* Parameters of a neuron, or of a homogeneous 
    neuron population. */

typedef struct nemasim_bundle_parms_t
  { double K;     /* Average number of synapses per source neuron. */
    double W_dev; /* Deviation of total rested synapse gain. */
    double W_avg; /* Average total rested synapse gain. */
  } nemasim_bundle_parms_t;
  /* Parameters of synapses from a homogeneous population {A} 
    of neurons to another homogeneous population {B}
    (possibly the same as {A}).  If the population {B}
    has {N} neurons, each neuron of {A} has a synapse
    to each neuron of {B} with independent probability 
    {K/N}.  The strength of that synapse, in its fully
    rested state, is a random variable  with log-normal 
    distribution, mean {W_avg/K}, and deviation {W_dev/K}. */

typedef struct nemasim_cohort_state_t
  { double V_avg; /* Average potential. */
    double V_dev; /* Deviation of potential.
    double G;     /* Input gain. */
    double H;     /* Output gain. */
    double eta;   /* Fraction of population that is in this cohort. */
    double rho;   /* Fraction of this cohort that will fire next. */
  } nemasim_cohort_state_t;
  /* State of a cohort of neurons with the same
    firing age {k} at some discrete time {t}.  
    The neurons have average potential {V_avg}
    with deviation {V_dev}, input gain {G}, output gain {H}, 
    and comprise the fraction {eta} of the population. 
    A fraction {rho} of them will fire during the 
    next time step.  */

typedef struct nemasim_pop_state_t
  { int32_t np;  /* Number of individual cohorts retained. */
    nemasim_cohort_state_t *cs; /* State of neurons in each cohort. */
  } nemansim_pop_state_t;
  /* Stochastic state of a homogeneous population of neurons.
    The neurons with firing age {k} have state {cs[k]},
    for {k} in {0..np-1}.  All neurons with age {np} or greater
    are lumped together and assumed to have state
    {cs[np].}  */

void nemasim_cohort_set
  ( nemasim_cohort_state_t *cs, 
    double V_avg,
    double V_dev,
    double G,
    double H,
    double eta,
    nemasim_parms_t *parm
  )
  /* Sets the state parameters of a cohort
    {cs} to potential {V_avg ± V_dev}, input gain factor {G},
    output gain factor {H}, and relative fraction
    {eta}.  The firing fraction {rho} is set 
    based on the firing function {parm->Phi}
    and the potential {V_avg,V_dev}. */
  {
    cs->V_avg = V_avg;
    cs->V_dev = V_dev;
    cs->G = G;
    cs->H = H;
    cs->eta = eta;
    cs->rho = nemasim_firing_prob(parm->Phi, V_avg, V_dev);
  }

void nemasim_cohort_evolve
  ( nemasim_cohort_state_t *cso, 
    double DV_avg,
    double DV_dev,
    double I_avg,
    double I_dev,
    nemasim_neuron_parms_t *parms, 
    nemasim_cohort_state_t *csn_fire,
    nemasim_cohort_state_t *csn_fail
  )
  /* Simulates the evolution of neurons in a specific cohort 
    {S[k,t]} of a homogeneous population, from some discrete time 
    {t} to the next time {t+1}.  The cohort is assumed to have some age 
    {k} and the stochastic state {cso} at time {t}.
    A fraction {cso->rho} of those neurons is assumed to fire, 
    and they will join the cohort {S[0,t+1]} of age zero at time
    {t+1}.  The sate of those neurons that fire is returned in {csn_fire}.
    The neurons that fail to fire
    will become the cohort {S[k+1,t+1]}, with age {k+1} at time {t+1},
    whose state is stored by the procedure into {csn_fail}.

    It is assumed that the total input that each neuron {i} 
    in {S[k,t]} received from all its input synapses betwen time {t}
    and time {t+1} is a normal random variable
    with mean {DV_avg} and deviation {DV_dev}.  This potential
    increment includes the pulses {X[j,t]} of each input neuron {j}
    in the network, modulated by its output gain factor {H[j,t]} and
    by the rested synaptic gain {wfix[j,i]}; but not
    yet modulated by the input gain {G[i,t]} of neuron {i}.

    In addition to synaptic inputs, each neuron {i} in the 
    cohort is also assumed to have received, between times
    {t} and {t+1}, an external input signal {I[i,t]},
    which is a normally distributed variable with mean
    {I_avg} and deviation {I_dev}.

    ??? Should the input signal be modulated by the input gain
    factor {G[i,t]} too? ???  */
  {
    /* Grab some variables for convenience: */
    double V_avg = cso->V_avg; /* Mean potential of cohort. */
    double V_dev = cso->V_dev; /* Potential deviation in cohort. */
    double G = cso->G;         /* Input gain factor of cohort. */
    double H - cso->H;         /* Output gain factor of cohort. */
    double eta = cso->eta;     /* Fraction of population in cohort. */
    double rho = cso->rho;     /* Fraction of cohort that is about to fire. */

    double mu_V = parm->mu_V;  /* Potential recovery factor. */
    double mu_G = parm->mu_G;  /* Input gain recovery factor. */
    double mu_H = parm->mu_H;  /* Output gain recovery factor. */

    /* Store the new state of the firing neurons in {csn_fire}: */
    nemasim_cohort_set
      ( csn_fire,
        parms->V_R, 0.0, parms->G_R, parms->H_R, eta*rho,
        parm
      );
   
    /* Compute the average potential of the non-firing neurons: */
    double V_avg_new = mu_V*V_avg + I_avg + DV_avg*G;
    /* Compute the deviation of the potential of those neurons: */
    /* ??? Can we ignore correlations here? ??? */
    double V_dev_new = hypot(hypot(mu_V*V_dev, I_dev), DV_dev*G);
    /* Compute the new input and  output gain factors of those neurons: */
    double G_new = 1.0 - mu_G*(1.0 - G);
    double H_new = 1.0 - mu_H*(1.0 - H);
    /* Store the new state of the non-firing neurons in {csn_fail}: */
    nemasim_cohort_set
      ( csn_fail,
        V_new_avg, V_new_dev, G_new, H_new, eta*(1.0 - rho),
        parm
      );
  }

void nemasim_pop_evolve
  ( nemasim_pop_state_t *pso,
    double DV_avg,
    double DV_dev,
    double I_avg,
    double I_dev,
    nemasim_neuron_parms_t *parm,
    nemasim_pop_state_t *psn
  )
  /* Simulates the evolution of a homogeneous 
    population of neurons with stochastic state {pso}
    and neuron parameters {parm}
    from some discrete time {t} to time {t+1}.

    It is assumed that the total input that each neuron {i} 
    in the population received from all its input synapses betwen time {t}
    and time {t+1} is a normal random variable
    with mean {DV_avg} and deviation {DV_dev}.  This potential
    increment includes the pulses {X[j,t]} of each input neuron {j}
    in the network, modulated by its output gain factor {H[j,t]} and
    by the rested synaptic gain {wfix[j,i]}; but not
    yet modulated by the input gain {G[i,t]} of neuron {i}.

    In addition to synaptic inputs, each neuron {i} in the 
    population is also assumed to have received, between times
    {t} and {t+1}, an external input signal {I[i,t]},
    which is a normally distributed variable with mean
    {I_avg} and deviation {I_dev}.

    ??? Should the input signal be modulated by the input gain
    factor {G[i,t]} too? ???  */
  {
    /* The state of the new cohort of age zero: */
    nemasim_cohort_state_t csn_fire;
    nemasim_cohort_clear(&csn_fire);
    /* The state of the new cohort of age {np}: */
    nemasim_cohort_state_t csn_lump;

    int k;
    for (k = np; k > 0; k--)
      { /* Evolve cohort {S[k,t]} to time {t+1}: */
        nemasim_cohort_state_t *cso = &(ps->cs[k]);
        nemasim_cohort_state_t csn_k_fire;
        nemasim_cohort_state_t csn_k_fail;
        nemasim_cohort_evolve
	  ( cso, DV_avg, DV_dev, I_avg, I_dev,
	    parms, 
	    &csn_k_fire, &csn_k_fail
	  );

        /* Merge the part of the cohort that fired into {csn_fire}: */
        nemasim_cohort_merge(&csn_k_fire, &csn_fire);

        /* Store or merge the part of the cohort that did not fire: */
        if (k == np)
          { /* Save the new cohort {np+1} to merge later: */
            csn_lump = csn_k_fail; 
	  }
        else 
          { nemasim_cohort_state_t *csn = &(ps->cs[k+1]); 
            (*csn) = csn_k_fail;
            if (k == np-1) 
               { nemasim_cohort_merge(&csn_lump, csn); }
	  }
      }
   
   }

void nemasim_network_evolve
  (
  )
  /* Simulates the evolution of a 
    network with one or more homogeneous 
    populations of neurons with stochastic states {net.pso[k]}
    and neuron parameters {parm} from some discrete time {t} 
    to time {t+1}.

    In addition to synaptic inputs, each neuron {i} in the 
    population {k} is also assumed to have received, between times
    {t} and {t+1}, an external input signal {I[i,t]},
    which is a normally distributed variable with mean
    and deviation given by {net.exin(k,t,&I_avg,I_dev)}

    ??? Should the input signal be modulated by the input gain
    factor {G[i,t]} too? ???  */
  {

    /* Compute the total input that each neuron {i} 
    in the population received from all its input synapses betwen time {t}
    and time {t+1} is a normal random variable
    with mean {DV_avg} and deviation {DV_dev}.  This potential
    increment includes the pulses {X[j,t]} of each input neuron {j}
    in the network, modulated by its output gain factor {H[j,t]} and
    by the rested synaptic gain {wfix[j,i]}; but not
    yet modulated by the input gain {G[i,t]} of neuron {i}. */


    /*  to the state of the cohort of 
    neurons with age 0 (those that just fired).
    Namely, the fields {cs.V,cs.G,cs.H} are set
    to the reset values {V_R,G_R,H_R}, respectively.
    The fraction {cs.eta} is set to {rho}. */
  }
