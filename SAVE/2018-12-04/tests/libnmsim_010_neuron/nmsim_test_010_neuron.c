#define PROG_NAME "nmsim_test_010_neuron"
#define PROG_DESC "basic tests of {limnmism} neuron-level procedures"
#define PROG_VERS "1.0"

/* Last edited on 2017-07-25 00:24:39 by stolfilocal */ 

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
#include <math.h>

#include <bool.h>
#include <argparser.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <affirm.h>

#include <nmsim_firing_func.h>
#include <nmsim_neuron_parms.h>
#include <nmsim_neuron_state.h>
#include <nmsim_neuron_net.h>

void nmsim_test_discrete_net(nmsim_neuron_count_t N, uint64_t T);
  /* Tests a neuron net with {N} neurons for {T} steps. */

nmsim_neuron_net_t *nmsim_test_discrete_net_new(nmsim_neuron_count_t N);
  /* Creates a neuron net with {N} neurons with uniform parameters
    and suitable random synapses. */
    
nmsim_neuron_parms_t *nmsim_test_discrete_net_parms_pick(void);
  /* Generates a set of neuron parameters suitable for testing. */

void nmsim_test_net_init(nmsim_neuron_net_t *net);
  /* Initializes the neurons of {net} with random states. */

void nmsim_test_net_simulate(nmsim_neuron_net_t *net, uint64_t T);
  /* Simulates the evolution of the network {net},
    from its current state, through {T} time steps. */

void nmsim_test_debug_firings(char *lab, nmsim_neuron_net_t *net, bool_t X[]);
void nmsim_test_debug_potentials(char *lab, nmsim_neuron_net_t *net);
void nmsim_test_debug_ages(char *lab, nmsim_neuron_net_t *net);
void nmsim_test_debug_potential_increments(char *lab, nmsim_neuron_net_t *net, double dV[]);
  /* Print state and state change variables to stderr, in one row. */

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
    
void nmsim_test_discrete_net(nmsim_neuron_count_t N, uint64_t T)
  {
    fprintf(stderr, "--- testing discrete network ---\n");
    fprintf(stderr, "network size {N} = %lu neurons\n", N);
    fprintf(stderr, "simulation length {T} = %lu steps\n", T);
    
    /* Creates a random network: */
    nmsim_neuron_net_t *net = nmsim_test_discrete_net_new(N);
    nmsim_test_net_init(net);
    nmsim_test_net_simulate(net, T);
  }
  
nmsim_neuron_net_t *nmsim_test_discrete_net_new(nmsim_neuron_count_t N)
  { 
    /* Create the neurons: */
    nmsim_neuron_net_t *net = nmsim_neuron_net_new(N);
    /* Set the neuron parameters: */
    nmsim_neuron_parms_t *parms = nmsim_test_discrete_net_parms_pick();
    for (nmsim_neuron_ix_t i = 0; i < net->N; i++) { net->parms[i] = parms; }
    
    /* Create random synapses: */
    double avg_deg = sqrt((double)N); /* Average in- and out-degree. */
    nmsim_synapse_count_t M = (nmsim_synapse_count_t)((double)N * floor(avg_deg));
    double wt_avg = 1.0/avg_deg;
    double wt_dev = 0.1*wt_avg;
    double pr_neg = 0.10; /* Probability of inhibitory synapse. */
    nmsim_neuron_net_create_synapses(net, M, wt_avg, wt_dev, pr_neg);
    return net;
  }
  
nmsim_neuron_parms_t *nmsim_test_discrete_net_parms_pick(void)
  { double V_B = -40.0;   
    double V_R = -30.0;   
    double c_B = 0.90;   
    double M_R = 1.00;   
    double M_mu = 1.00;  
    double G_R = 1.00;   
    double G_mu = 1.00;  
    double H_R = 1.00;   
    double H_mu = 1.00;  
    
    /* Firing function parameters: */
    double V_M = -20.0;
    double D_M = 2.0;
    nmsim_firing_func_t *Phi = nmsim_firing_func_gauss_new(V_M, D_M);
    nmsim_neuron_parms_t *parms = 
      nmsim_neuron_parms_new(V_B, V_R, c_B, M_R, M_mu, G_R, G_mu, H_R, H_mu, Phi);
    return parms;
  }

void nmsim_test_net_init(nmsim_neuron_net_t *net)
  {
    for (nmsim_neuron_ix_t i = 0; i < net->N; i++) 
      { nmsim_neuron_state_t *state = &(net->state[i]); 
        nmsim_neuron_parms_t *parms = net->parms[i];
        double rho = 0.005; /* For now, assume firing every 100 steps. */
        nmsim_neuron_state_throw(state, parms, rho);
      }
  }

void nmsim_test_net_simulate(nmsim_neuron_net_t *net, uint64_t T)
  {
    nmsim_neuron_count_t N = net->N;
    bool_t debug = FALSE;
    bool_t *X = notnull(malloc(N*sizeof(bool_t)), "no mem"); /* Firing indicators. */ 
    double *dV = notnull(malloc(N*sizeof(double)), "no mem"); /* INput potential increments. */ 
    for (uint64_t t = 0; t < T; t++)
      { 
        if (debug) { fprintf(stderr, "iteration %lu -> %lu\n", t, t+1); }

        if (debug) { nmsim_test_debug_potentials("V0", net); }
        /* Fisrt pass: compute the firing indicators: */
        for (nmsim_neuron_ix_t i = 0; i < net->N; i++) 
          { nmsim_neuron_parms_t *parms = net->parms[i];
            double pr;
            nmsim_firing_func_eval(parms->Phi, net->state[i].V, &pr, NULL);
            X[i] = (drandom() < pr);
          }
        if (debug) { nmsim_test_debug_firings("X  = ", net, X); }

        /* Reset or decay the previous potentials: */
        for (nmsim_neuron_ix_t i = 0; i < net->N; i++) 
          { nmsim_neuron_parms_t *parms = net->parms[i];
            nmsim_neuron_state_t *state = &(net->state[i]);
            if (X[i]) 
              { state->V = parms->V_R; }
            else
              { double c = nmsim_neuron_state_leakage(parms, state->age);
                state->V = c * state->V;
              }
          }
        if (debug) { nmsim_test_debug_potentials("V1", net); }

        /* Integrate the inputs: */
        for (nmsim_neuron_ix_t i = 0; i < net->N; i++) 
          { nmsim_neuron_state_t *state = &(net->state[i]);
            dV[i] = nmsim_neuron_net_tot_input(net, X, i);
            state->V += dV[i];
         }
        if (debug) { nmsim_test_debug_potential_increments("dV", net, dV); }
        if (debug) { nmsim_test_debug_potentials("V2", net); }

        /* Update the firing ages: */
        if (debug) { fprintf(stderr, "V1 = "); }
        for (nmsim_neuron_ix_t i = 0; i < net->N; i++) 
          { nmsim_neuron_state_t *state = &(net->state[i]);
            state->age = (X[i] ? 0 : state->age + 1);  
            if (debug) { fprintf(stderr, " %6.2f", state->V); }         
          }
        if (debug) { fprintf(stderr, "\n"); }
      }
  }

void nmsim_test_debug_firings(char *lab, nmsim_neuron_net_t *net, bool_t X[])
  { fprintf(stderr, "%s = ", lab);        
    for (nmsim_neuron_ix_t i = 0; i < net->N; i++) 
      { fprintf(stderr, " %6d", (int)(X[i])); }
    fprintf(stderr, "\n");
  }
          
void nmsim_test_debug_potentials(char *lab, nmsim_neuron_net_t *net)
  { fprintf(stderr, "%s = ", lab);        
    for (nmsim_neuron_ix_t i = 0; i < net->N; i++) 
      { nmsim_neuron_state_t *state = &(net->state[i]);
        fprintf(stderr, " %6.2f", state->V);
      }
    fprintf(stderr, "\n");
  }
           
void nmsim_test_debug_ages(char *lab, nmsim_neuron_net_t *net)
  { fprintf(stderr, "%s = ", lab);        
    for (nmsim_neuron_ix_t i = 0; i < net->N; i++) 
      { nmsim_neuron_state_t *state = &(net->state[i]);
        fprintf(stderr, " %6lu", state->age);
      }
    fprintf(stderr, "\n");
  }
           
void nmsim_test_debug_potential_increments(char *lab, nmsim_neuron_net_t *net, double dV[])
  { fprintf(stderr, "%s = ", lab);        
    for (nmsim_neuron_ix_t i = 0; i < net->N; i++) 
      { fprintf(stderr, " %6.2f", dV[i]); }
    fprintf(stderr, "\n");
  }
         
