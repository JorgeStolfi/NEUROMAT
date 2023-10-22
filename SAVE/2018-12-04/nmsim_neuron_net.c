/* See {nmsim_neuron_net.h} */
/* Last edited on 2017-07-25 00:24:23 by stolfilocal */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <jsrandom.h>

#include <nmsim_neuron_parms.h>
#include <nmsim_neuron_state.h>

#include <nmsim_neuron_net.h>

nmsim_neuron_net_t *nmsim_neuron_net_new(nmsim_neuron_count_t N)
  { 
    demand(N <= nmsim_neuron_net_MAX_NEURONS, "too many neurons");
    
    nmsim_neuron_net_t *net = notnull(malloc(sizeof(nmsim_neuron_net_t)), "no mem");
    net->N = N;
    net->state = notnull(malloc(N*sizeof(nmsim_neuron_state_t)), "no mem");
    net->parms = notnull(malloc(N*sizeof(nmsim_neuron_parms_t *)), "no mem");
    net->deg_in = notnull(malloc(N*sizeof(nmsim_neuron_count_t)), "no mem");
    net->src_in = NULL;
    net->deg_ot = notnull(malloc(N*sizeof(nmsim_neuron_count_t)), "no mem");
    net->dst_ot = NULL;
    for (nmsim_neuron_ix_t i = 0; i < N; i ++)
      { net->deg_in[i] = 0;
        net->deg_ot[i] = 0;
        net->parms[i] = NULL;
        net->state[i] = (nmsim_neuron_state_t){ .V = NAN, .age = (~(uint64_t)0) };
      }
    return net;
  }

void nmsim_neuron_net_create_synapses
  ( nmsim_neuron_net_t *net,
    nmsim_synapse_count_t M,
    double wt_avg,
    double wt_dev,
    double pr_neg
  )
  {
    nmsim_neuron_count_t N = net->N;
    /* Clear degree counts, pick inhibitors, make sure that synapse lists are {NULL}: */
    bool_t *inhib = notnull(malloc(N*sizeof(bool_t)), "no mem"); /* Inhibitory flag. */
    for (nmsim_neuron_ix_t i = 0; i < N; i++)
      { net->deg_in[i] = 0;
        net->deg_ot[i] = 0;
        if (net->src_in[i] != NULL) { free(net->src_in[i]); net->src_in[i] = NULL; }
        if (net->dst_ot[i] != NULL) { free(net->dst_ot[i]); net->dst_ot[i] = NULL; }
        inhib[i] = (drandom() <= pr_neg);
      }
    /* Generate the synapses as pairs {pre[k],pos[k]} for {k} in {0..M-1}: */
    /* Also accumulate degrees: */
    nmsim_neuron_ix_t *pre = notnull(malloc(M*sizeof(nmsim_neuron_ix_t)), "no mem"); /* Pre-synaptic neurons. */
    nmsim_neuron_ix_t *pos = notnull(malloc(M*sizeof(nmsim_neuron_ix_t)), "no mem"); /* Post-synaptic neurons. */
    for (nmsim_synapse_ix_t k = 0; k < M; k++)
      { nmsim_neuron_ix_t j = uint64_abrandom(0, N); /* Pre-synaptic (source) neuron. */
        nmsim_neuron_ix_t i = uint64_abrandom(0, N); /* st-synaptic (dest) neuron. */
        pre[k] = j; net->deg_ot[j]++;
        pos[k] = i; net->deg_in[i]++;
      }
    /* Second pass: Allocate input and output synapse tables, reset degrees: */
    for (nmsim_neuron_ix_t i = 0; i < net->N; i++) 
      { /* Allocate input source and weight table: */
        nmsim_neuron_count_t din = net->deg_in[i];
        assert(din <= M);
        nmsim_neuron_ix_t *src = notnull(malloc(din*sizeof(nmsim_neuron_ix_t)), "no mem");
        double *wt = notnull(malloc(din*sizeof(nmsim_neuron_ix_t)), "no mem");
        for(nmsim_synapse_ix_t k = 0; k < din; k++) { src[k] = i; wt[k] = 0.0; /* For now. */ }
        net->src_in[i] = src;
        net->wt_in[i] = wt;
        net->deg_in[i] = 0;

        /* Allocate output destination table: */
        nmsim_neuron_count_t dot = net->deg_ot[i];
        assert(dot <= M);
        nmsim_neuron_ix_t *dst = notnull(malloc(dot*sizeof(nmsim_neuron_ix_t)), "no mem");
        for(nmsim_synapse_ix_t k = 0; k < dot; k++) { dst[k] = i; /* For now. */ }
        net->dst_ot[i] = dst;
        net->deg_ot[i] = 0;
      }

    /* Third pass:  Collect synapses, set weights, increment degrees. */
    double wt_avg_log = log(wt_avg);
    double wt_dev_log = log((wt_avg + wt_dev)/wt_avg);
    for (nmsim_synapse_ix_t k = 0; k < M; k++)
      { /* Grab neuron indices: */
        nmsim_neuron_ix_t j = pre[k]; /* Pre-synaptic (source) neuron. */
        nmsim_neuron_ix_t i = pos[k]; /* Post-synaptic (dest) neuron. */
        
        /* Compute the synapse's rest weight {wtk}, store as the input weight of {j}: */
        double wtk = exp(wt_avg_log + wt_dev_log*dgaussrand());
        if (inhib[i]) { wtk = - wtk; }

        /* Store indices and weight in the synapse tables: */
        nmsim_neuron_ix_t ki = net->deg_in[i]; 
        net->src_in[i][ki] = j;
        net->wt_in[i][ki] = wtk;
        net->deg_in[i]++;
        
        nmsim_neuron_ix_t kj = net->deg_ot[j];
        net->dst_ot[j][kj] = i;
        net->deg_ot[j]++;
      }
      
    /* Free auxiliary storage: */
    free(inhib);
    free(pre);
    free(pos);
  }

double nmsim_neuron_net_tot_input(nmsim_neuron_net_t *net, bool_t X[], nmsim_neuron_ix_t i)
  { 
    /* Grab properties of neuron {i}: */
    nmsim_neuron_parms_t *parms_i = net->parms[i];
    nmsim_neuron_state_t *state_i = &(net->state[i]);
    nmsim_neuron_count_t deg_in_i = net->deg_in[i]; /* Number of input synapses. */
    nmsim_neuron_ix_t *src_i = net->src_in[i]; /* Indices of presynaptic neurons. */
    double *wt_i = net->wt_in[i];  /* Resting weights of input synapses. */
    double G_i = nmsim_neuron_state_input_gain(parms_i, state_i->age);
    
    /* Scan the input synapses of neuron {i}: */
    double dV = 0.0;
    for (nmsim_synapse_ix_t ki = 0; ki < deg_in_i; ki++) 
      { /* Get the index of {j} of the pre-synaptic neuron: */
        nmsim_neuron_ix_t j = src_i[ki];
        if (X[j])
          { /* Grab parameters of neuron {j}: */
            nmsim_neuron_parms_t *parms_j = net->parms[j];
            nmsim_neuron_state_t *state_j = &(net->state[j]);
            double H_j = nmsim_neuron_state_output_gain(parms_j, state_j->age);
            double wt_ji = wt_i[ki]; /* Resting weight of sunapse. */
            /* Accumulate the input from neuron {j}: */
            dV += G_i * H_j * wt_ji;
          }
      }
    return dV;
  }
