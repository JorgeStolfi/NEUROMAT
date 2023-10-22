#define PROG_NAME "nmsim_elem_net_simulate"
#define PROG_DESC "element-level simulation of a GL neuron network"
#define PROG_VERS "1.0"

/* Last edited on 2020-12-24 23:32:48 by jstolfi */ 

#define PROG_COPYRIGHT \
  "Copyright © 2019  State University of Campinas (UNICAMP)"
  
#define PROG_AUTH \
  "Created by J. Stolfi on 2019-03-21"
  
#define PROG_HIST
  
#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -netFile {netFile} \\\n" \
  "    [ -timeStep {timeStep} ] \\\n" \
  "    -nSteps {nSteps} \\\n" \
  "    [ -exInput {ineLo} {ineHi} tLo} {tHi} {Iavg} {Idev} ].. \\\n" \
  "    [ -trace {ineLo} {ineHi} tLo} {tHi} ].. \\\n" \
  "    [ -groupStats tLo} {tHi} ] \\\n" \
  "    -outPrefix {outPrefix} \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program reads an element-level description of" \
  " a Galves-Loecherbach neuron net and simulates its" \
  " behavior for a specified number of steps, with specified" \
  " external inputs. \n" \
  "\n" \
  "  The program writes several files \"{outPrefix}_elem_*_trace.txt\" with" \
  " traces of selected neurons in selected time" \
  " intervals (see the \"-trace\" option).  It also writes" \
  " a statstical summary of the network's synapses to the" \
  " file \"{outPrefix}_synapse_stats.txt\".\n" \
  "\n" \
  "OPTIONS\n" \
  "  -netFile {netFile} \n" \
  "    This mandatory command line argument specifies the" \
  " file name {netFile} containing the description of the" \
  " network.\n" \
  "\n" \
  "  -outPrefix {outPrefix} \n" \
  "    This mandatory argument specifies the common" \
  " prefix {outPrefix} for all output file names.\n" \
  "\n" \
  "\n" \
  "  -timeStep {timeStep} \n" \
  "    This optional argument specifies the simulation" \
  " time step length {timeStep} in milliseconds.  If omitted," \
  " assumes \"-timeStep 1\" (1 millisecond).\n" \
  "\n" \
  "  -nSteps {nSteps} \n" \
  "    This mandatory argument specifies the" \
  " number {nSteps} of time steps to simulate.  The" \
  " discrete times will range from 0 to {nSteps}.  (However, the" \
  " simulation parameters {X,I,J} will not be defined for the last time {nSteps}.)\n" \
  "\n" \
  "  -exInput {ineLo} {ineHi} {tLo} {tHi} {Iavg} {Idev}\n" \
  "    This argument specifies the next-step external input" \
  " increment {I[k][t]} for the neurons with indices {k} in" \
  " the range {ineLo .. ineHi}, for the discrete" \
  " times {t} in the range {tLo .. tHi}. It can be" \
  " repeated as many times as desired.\n" \
  "\n" \
  "    If the same" \
  " pair {k,t} is specified in two or more of these" \
  " specs, the corresponding values are added together.  For any" \
  " neurons and times not specified via these arguments, the" \
  " external input {I[k][t]} will be zero.\n" \
  "\n" \
  "    The unit of {Iavg} is millivot per millisecond, and the unit of {Idev} is millivolt" \
  " per square root of millisecond. The inputs" \
  " {I[k][t]} will be generated independently of each" \
  " other and of any other inputs or parameters, by sampling a Gaussian" \
  " distribution with mean {Iavg*timeStep}, and deviation {Idev*sqrt(timeStep)}, where" \
  " {timeStep} is expressed in milliseconds.\n" \
  "\n" \
  "  -trace {ineLo} {ineHi} {tLo} {tHi}\n" \
  "    This optional argument specifies that the program should write" \
  " out the average of states and activities of neurons with indices" \
  " in {ineLo..ineHi} for" \
  " each discrete time {t} in {tLo..tHi}.  In particular, if" \
  " {ineLo=ineHi}, the state and firings of that single neuron" \
  " is written out. This argument can specified" \
  " more than once, to trace different neurons or neuron sets over different time" \
  " intervals; but the neuron index ranges had better be disjoint.   The data" \
  " for each use of this argument is written to a separate file" \
  " called \"{outPrefix}_ne{ineLo}--{ineHi}.txt\" where the indices {ineLo} and {ineHi}" \
  " are formatted as 10 decimal digits, zero-padded.\n" \
  "\n" \
  "  -groupStats {tLo} {tHi}\n" \
  "    This optional argument specifies that the program should compute" \
  " statistics of the state variables (external inputs, potential, firing" \
  " rate, etc.) for all neurons in each neuron group and all discrete" \
  " times {t} in {tLo..tHi}.  This argument can specified" \
  " only once.  " nmsim_elem_net_sim_group_stats_write_INFO "\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  nmsim_group_net_simulate(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2013-06-01 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2019-03-21 J. Stolfi: created.\n" \
  "  2020-12-06 J. Stolfi: scaled external input parameters with time step.\n" \
  "  2020-12-14 J. Stolfi: added synapse parameter statistics by group.\n" \
  "  2020-12-15 J. Stolfi: added \"-groupStats\" - state variable statistics by group.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " PROG_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define stringify(x) strngf(x)
#define strngf(x) #x

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
#include <jsmath.h>
#include <jstime.h>
#include <affirm.h>

#include <nmsim_basic.h>
#include <nmsim_class_net.h>
#include <nmsim_group_net.h>
#include <nmsim_elem_net.h>
#include <nmsim_elem_net_group_stats.h>
#include <nmsim_elem_net_trace.h>

#include <nmsim_elem_net_sim.h>

typedef struct nens_elem_neuron_time_range_t
  { nmsim_elem_neuron_ix_t ineLo; /* Index of first neuron in set. */
    nmsim_elem_neuron_ix_t ineHi; /* Index of last neuron in set. */
    nmsim_time_t tLo;           /* First time in range. */
    nmsim_time_t tHi;           /* Last time in range. */
  } nens_elem_neuron_time_range_t;
  /* Specifies a range of neuron indices {ineLo..ineHi} 
    and a range of times {tLo..tHi}. */

typedef struct nens_input_spec_t
  { nens_elem_neuron_time_range_t kt;    /* Range of neurons and times. */
    double I_avg;                        /* Mean input value (mV). */
    double I_dev;                        /* Deviation of input values (mV). */
  } nens_input_spec_t;
  /* Description of external inputs for a range of neurons
    during a range or times.  Each input {I[k][t]} will be 
    a random variable independently drawn from a Gaussian 
    distribution with mean {I_avg} and deviation {I_dev}. 
    Note that these values are assumed to be adjusted 
    for the time step. */
    
vec_typedef(nens_input_spec_vec_t,nens_input_spec_vec,nens_input_spec_t);

typedef struct nens_trace_spec_t
  { nens_elem_neuron_time_range_t kt;    /* Range of neurons and times. */
  } nens_trace_spec_t;
  /* Specification of a range of neurons whose evolution
    should be monitored during a specified time range. */
    
vec_typedef(nens_trace_spec_vec_t,nens_trace_spec_vec,nens_trace_spec_t);
    
typedef struct nens_gstats_spec_t
  { nmsim_time_t tLo;    /* First time to gather stats. */
    nmsim_time_t tHi;    /* Last time to gather stats. */
  } nens_gstats_spec_t;
  /* Specification of a range of times for collecting per-grup statistics of neuron
    state and activity variables. */
    
typedef struct nens_options_t
  { char *netFile;                    /* Name of file with the network's description. */
    double timeStep;                  /* Step length to use in simulation (ms). */
    nmsim_step_count_t nSteps;        /* Number of steps to simulate. */
    nens_input_spec_vec_t exInput;    /* Specs of external neuron inputs (scaled for {timeStep)}. */
    nens_trace_spec_vec_t trace;      /* Neurons and times to trace. */
    nens_gstats_spec_t groupStats;    /* Specs for collecting neuron stats by group. */
    char *outPrefix;                  /* Prefix for output file names. */
  } nens_options_t;
  /* Arguments from command line. */
 
nmsim_elem_net_t *nens_read_network(char *fname, double  timeStep);
  /* Reads a descriptin of an element-level GL neuron net 
    from file {rd}. */

void nens_write_synapse_summary(char *prefix, nmsim_elem_net_t *enet);
  /* Writes to file "{prefix}_synapse_stats.txt" a summary of the synapse 
    parameters of {enet} by synapse group. */
    
nmsim_elem_net_trace_t *nens_make_trace
  ( nmsim_elem_net_t *enet, 
    nmsim_time_t sim_tLo,
    nmsim_time_t sim_tHi,
    nens_trace_spec_vec_t trspec
  );
  /* Creates a trace data structure for the neuron net {enet},
    that will store the full states of neurons and times 
    specified in {trspec}, clipped to the time range {sim_tLo..sim_tHi}. */
    
nmsim_elem_net_sim_group_stats_t *nens_make_gstats
  ( nmsim_elem_net_t *enet, 
    nmsim_time_t sim_tLo, 
    nmsim_time_t sim_tHi,
    nens_gstats_spec_t *gsspecs
  );
  /* Creates a gata structure to store the per-group neuron state and activiviy
    statistics for the network {enet},
    covering the time range specified in {gsspecs}
    clipped to the simulation range {sim_tLo..sim_tHi}. */

typedef void nens_trace_spec_visit_proc_t 
  ( nmsim_elem_neuron_ix_t ineLo_r, 
    nmsim_elem_neuron_ix_t ineHi_r, 
    nmsim_time_t tLo_r, 
    nmsim_time_t tHi_r
  );
  /* Type of a procedure that processes a range of neuron indices
    {ineLo_r..ineHi_r} and a range of times {tLo_r..tHi_r},
    called by {nens_scan_trace_specs} below.  The argument
    ranges will be non-empty and contained in the global 
    index and time intervals. */
    
void nens_scan_trace_specs
  ( nens_trace_spec_vec_t trspec,
    nmsim_elem_neuron_count_t nne,
    nmsim_time_t sim_tlo,
    nmsim_time_t sim_tHi,
    bool_t verbose,
    nens_trace_spec_visit_proc_t *visit
  );
  /* Scans the trace specifications {trspec}.  For each entry
    whose neuron index interval intersects the global index interval {0..nne-1} and 
    whose time interval intersects the global time interval {sim_tLo..sim_tHi}, calls {visit}
    with the intersections of those intervals.
    
    If {verbose} is true, prints warnings to {stderr}
    when the {trspec} entry intervals are not contained in the
    global intervals. */

void nens_set_external_inputs
  ( nmsim_elem_neuron_count_t nne, 
    double I[], 
    nmsim_time_t t,
    nens_input_spec_vec_t exInput
  );
  /* Defines the external inputs {I[k]} for the time step from {t} to {t+1},
    for all neurons {k} in {0..nne-1}, according to the specs in {exInput}.  Namely,
    for entry {exr} of {exInput} such that {t} is in the 
    time range {exr.tLo .. exr.tHi}, and for every neuron index 
    {k} in {exr.ineLo .. exr.ineHi}, sets {I[k]} to a log-normal random variable
    with the mean and deviation specified in {exInput}. If the same pair {k,t} is specified
    multiple times in {exInput}, the corresponding values are added together. 
    For any neurons and times not specified in {exInput}, the value {I[k]} is set to zero.  */

void nens_simulate
  ( nmsim_elem_net_t *enet, 
    nens_input_spec_vec_t exInput,
    nmsim_time_t tLo,
    nmsim_time_t tHi,
    nmsim_elem_net_trace_t *etrace,
    nmsim_elem_net_sim_group_stats_t *gstats
  );
  /* Simulates the evolution of the network {enet} from time {t=tLo}
    to time {t=tHi} -- that is, for {tHi-tLo} time steps.  
    
    Also stores
    into {etrace} the full states of some neurons in some time ranges,
    as specified in {etrace}. 
    
    Also stores into {gstats} the per-group
    neuron state and activity statistics in the time range 
    specfied in {gstats}.
    
    Note that the states for {t=tHi} are only partially defined. */
  
nmsim_elem_net_t *nens_read_network(char *fname, double timeStep);
  /* Reads a network description from {fname}.  The {timeStep}
    is used to convert characteristc decay times to decay factors. */

void nens_write_traces(char *prefix, nmsim_elem_net_trace_t *etrace);
  /* Writes the elem-level traces in {etrace} with names "{prefix}_elem_n{NN}.txt"
    where {NN} is the neuron index. */

void nens_write_group_stats(char *prefix, nmsim_elem_net_sim_group_stats_t *gstats);
  /* Writes the per-group neuron state statistics collected {gstats}. 
    See {nmsim_elem_net_sim_group_stats_write}. */

nens_options_t *nens_parse_options(int argc, char **argv);
  /* Parses the command line arguments. */

nens_elem_neuron_time_range_t parse_elem_neuron_time_range(argparser_t *pp);
  /* Parses four arguments {ILO} {IHI} {tLo} {tHi} 
    as a range of neuron indices {ILO..IHI} and a range
    of times {tLo..tHi}. */

int main(int argc, char **argv);

/* IMPLEMENTATIONS: */

int main(int argc, char **argv)
  { 
    nens_options_t *o = nens_parse_options(argc, argv);
    nmsim_elem_net_t *enet = nens_read_network(o->netFile, o->timeStep);
    nens_write_synapse_summary(o->outPrefix, enet);
    nmsim_time_t sim_tLo = 0;                 /* Time at start of first step of simulation. */
    nmsim_time_t sim_tHi = sim_tLo + o->nSteps; /* Time at end of last step of simulation. */
    nmsim_elem_net_trace_t *etrace = nens_make_trace(enet, sim_tLo, sim_tHi, o->trace);
    nmsim_elem_net_sim_group_stats_t *gstats = nens_make_gstats(enet, sim_tLo, sim_tHi, &(o->groupStats));
    nens_simulate(enet, o->exInput, sim_tLo, sim_tHi, etrace, gstats);
    nens_write_traces(o->outPrefix, etrace);
    nens_write_group_stats(o->outPrefix, gstats);
    nmsim_elem_net_trace_free(etrace);
    return 0;
  }
  
void nens_simulate
  ( nmsim_elem_net_t *enet, 
    nens_input_spec_vec_t exInput,
    nmsim_time_t tLo,
    nmsim_time_t tHi,
    nmsim_elem_net_trace_t *etrace,
    nmsim_elem_net_sim_group_stats_t *gstats
  )
  {
    nmsim_group_net_t *gnet = enet->gnet; /* Description of neuron and synapse groups. */
    nmsim_class_net_t *cnet = gnet->cnet; /* Description of neuron and synapse classes. */
    
    nmsim_elem_neuron_count_t nne = enet->nne; /* Number of neurons. */
    
    /* Allocate work arrays: */
    fprintf(stderr, "Allocating work arrays for %d neurons...\n", nne);
    double *V = notnull(malloc(nne*sizeof(double)), "no mem");
    nmsim_step_count_t *age = notnull(malloc(nne*sizeof(nmsim_step_count_t)), "no mem");
    bool_t *X = notnull(malloc(nne*sizeof(bool_t)), "no mem");
    double *M = notnull(malloc(nne*sizeof(double)), "no mem");
    double *H = notnull(malloc(nne*sizeof(double)), "no mem");
    double *I = notnull(malloc(nne*sizeof(double)), "no mem");
    double *J = notnull(malloc(nne*sizeof(double)), "no mem");
    
    /* Initialize the neuron ages, potentials, and modulators: */
    fprintf(stderr, "initializing network state...\n");
    for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++) 
      { nmsim_group_neuron_ix_t ing = enet->neu[ine].ing; /* Neuron group index. */
        nmsim_class_neuron_ix_t inc = gnet->ngrp[ing].inc; /* Neuron class index. */
        nmsim_class_neuron_throw_state(cnet->nclass[inc], &(V[ine]), &(age[ine]));
      }
    nmsim_elem_net_sim_compute_modulators(enet, tLo, age, M, H);
    
    /* Initialize group-level statistics records {ngs[0..nng-1]}: */
    nmsim_elem_net_sim_group_stats_initialize(gstats);
    
    /* Simulate: */
    fprintf(stderr, "simulating for time in %ld..%ld ...\n", tLo, tHi);
    double rt0 = user_cpu_time_usec();
    for (nmsim_time_t t = tLo; t < tHi; t++)
      { if ((t - tLo) % 100 == 0) { fprintf(stderr, "."); }
        /* Define the external inputs for tHis time step: */
        nens_set_external_inputs(nne, I, t, exInput);
        /* Apply the evolution equations: */
        nmsim_elem_net_sim_step(enet, t, V, age, M, H, X, I, J, etrace, gstats);
      }
    /* Save the non-step data for the last time {t=tHi}: */
    if (etrace != NULL) { nmsim_elem_net_trace_set_V_age_M_H(etrace, tHi, nne, V, age, M, H); }  
    if (gstats != NULL) 
      { nmsim_elem_net_sim_group_stats_accumulate_V_age_M_H(gstats, tHi, nne, V, age, M, H); 
        nmsim_elem_net_sim_group_stats_finalize(gstats);
      }

    fprintf(stderr, "\n");
    double rt1 = user_cpu_time_usec();
    fprintf(stderr, "computer time for simulation = %.3f s...\n", (rt1-rt0)/1000000);
    
    /* Reclaim storage: */
    free(V);
    free(age);
    free(X);
    free(M); free(H);
    free(I); free(J);
  }
    
void nens_write_synapse_summary(char *prefix, nmsim_elem_net_t *enet)
  {
    char *fname = NULL;
    asprintf(&fname, "%s_synapse_stats.txt", prefix);
    FILE *wr = open_write(fname, TRUE);
    nmsim_elem_net_group_stats_print_all(wr, enet, "  ", "\n");
    fclose(wr);
    free(fname);
  }
    
void nens_scan_trace_specs
  ( nens_trace_spec_vec_t trspec,
    nmsim_elem_neuron_count_t nne,
    nmsim_time_t sim_tLo,
    nmsim_time_t sim_tHi,
    bool_t verbose,
    nens_trace_spec_visit_proc_t *visit
  )
  { 
  }

nmsim_elem_net_trace_t *nens_make_trace
  ( nmsim_elem_net_t *enet, 
    nmsim_time_t sim_tLo, 
    nmsim_time_t sim_tHi,
    nens_trace_spec_vec_t trspec
  )
  {
    bool_t verbose = FALSE;
    
    nmsim_elem_neuron_count_t nne = enet->nne; /* Num of neurons in network. */
        
    /* Count neuron sets to be traced: */
    
    nmsim_elem_neuron_count_t ntr = trspec.ne; /* Number of neurons to trace. */
    
    /* Create trace structure for the specified neuron seets, initially with empty time range: */
    fprintf(stderr, "creating trace structure for %d neuron sets...\n", ntr);
    nmsim_elem_net_trace_t *etrace = nmsim_elem_net_trace_new(ntr);
    for (int32_t ktr = 0; ktr < ntr; ktr++)
      { nens_trace_spec_t *spk = &(trspec.e[ktr]); /* Next neuron & time range spec. */
        /* Get and clip the specified neuron index range: */
        nmsim_elem_neuron_ix_t ineLo_clp, ineHi_clp;
        char *iname = (verbose ? "neuron index" : NULL);
        nmsim_int32_range_clip(spk->kt.ineLo, spk->kt.ineHi, 0, nne-1, &ineLo_clp, &ineHi_clp, iname); 
        /* Get and clip the time range: */
        nmsim_time_t tLo_clp, tHi_clp;
        char *tname = (verbose ? "time" : NULL);
        nmsim_int64_range_clip(spk->kt.tLo, spk->kt.tHi, sim_tLo, sim_tHi, &tLo_clp, &tHi_clp, tname); 
        /* Allocate the record if not empty: */
        bool_t empty = FALSE;
        if (ineLo_clp > ineHi_clp)
          { fprintf(stderr, "!! warning: neuron range %d .. %d", spk->kt.ineLo, spk->kt.ineHi);
            fprintf(stderr, " fully outside valid range\n"); 
            empty = TRUE;
          }
        if (tLo_clp > tHi_clp)
          { fprintf(stderr, "!! warning: time range %ld .. %ld", spk->kt.tLo, spk->kt.tHi);
            fprintf(stderr, " fully outside the simulated range\n"); 
            empty = TRUE;
          }
        if (! empty)
          { nmsim_elem_neuron_trace_t *trne = 
              nmsim_elem_neuron_trace_new(ineLo_clp, ineHi_clp, tLo_clp, tHi_clp);
            nmsim_elem_net_trace_set(etrace, ktr, trne);
          }
      }
    
    return etrace;
  }

nmsim_elem_net_sim_group_stats_t *nens_make_gstats
  ( nmsim_elem_net_t *enet, 
    nmsim_time_t sim_tLo, 
    nmsim_time_t sim_tHi,
    nens_gstats_spec_t *gsspecs
  )
  { nmsim_time_t ngs_tLo = (nmsim_time_t)imax(sim_tLo, gsspecs->tLo);
    nmsim_time_t ngs_tHi = (nmsim_time_t)imin(sim_tHi, gsspecs->tHi);
    nmsim_elem_net_sim_group_stats_t *gstats = nmsim_elem_net_sim_group_stats_new(enet->gnet, ngs_tLo, ngs_tHi);
    return gstats;
  }

void nens_set_external_inputs
  ( nmsim_elem_neuron_count_t nne, 
    double I[], 
    nmsim_time_t t,
    nens_input_spec_vec_t exInput
  )
  {
    bool_t debug = FALSE; /* (t == 400); */

    if (debug) { fprintf(stderr, "setting external inputs of neurons %d..%d\n", 0, nne-1); }

    /* Clear the inputs: */
    for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++) { I[ine] = 0.0; } 
    /* Get the external inputs for this time step: */
    for (int32_t r = 0; r < exInput.ne; r ++)
      { /* Grab the input spec: */
        nens_input_spec_t *exr = &(exInput.e[r]);
        if ((t >= exr->kt.tLo) && (t <= exr->kt.tHi))
          { nmsim_elem_neuron_ix_t ineLo = (nmsim_elem_neuron_ix_t)imax(0, exr->kt.ineLo);
            nmsim_elem_neuron_ix_t ineHi = (nmsim_elem_neuron_ix_t)imin(nne-1, exr->kt.ineHi);
            double Iavg = exr->I_avg;
            double Idev = exr->I_dev;
            if (debug) 
              { fprintf(stderr, "  adding %+8.4f +/- %8.4f to I[%d..%d]\n", Iavg, Idev, ineLo, ineHi); }
            for (nmsim_elem_neuron_ix_t ine = ineLo; ine <= ineHi; ine++) 
              { I[ine] += Iavg + Idev * dgaussrand(); }
          }
      }

    if (debug) 
      { for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++)
          { fprintf(stderr, "  external input of neuron %d = %10.4f\n", ine, I[ine]); }
      }
  }

nmsim_elem_net_t *nens_read_network(char *fname, double timeStep)
  { fprintf(stderr, "reading network from file %s\n", fname);
    FILE *rd = open_read(fname, TRUE);
    double rt0 = user_cpu_time_usec();
    nmsim_elem_net_t *enet = nmsim_elem_net_read(rd, timeStep);
    double rt1 = user_cpu_time_usec();
    fprintf(stderr, "computer time reading = %.3f s\n", (rt1-rt0)/1000000);
    fclose(rd);
    return enet;
  }

void nens_write_traces(char *prefix, nmsim_elem_net_trace_t *etrace)
  { char *epref = NULL;
    asprintf(&epref, "%s_elem", prefix);
    nmsim_elem_net_trace_write(epref, etrace);
    free(epref);
  }
    
void nens_write_group_stats(char *prefix, nmsim_elem_net_sim_group_stats_t *gstats)
  { nmsim_elem_net_sim_group_stats_write(prefix, gstats);
  }

vec_typeimpl(nens_input_spec_vec_t,nens_input_spec_vec,nens_input_spec_t);
    
vec_typeimpl(nens_trace_spec_vec_t,nens_trace_spec_vec,nens_trace_spec_t);
         
nens_options_t *nens_parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    /* Allocate the command line argument record: */
    nens_options_t *o = notnull(malloc(sizeof(nens_options_t)), "no mem");

    /* Parse keyword parameters: */

    argparser_get_keyword(pp, "-netFile");
    o->netFile = argparser_get_next_non_keyword(pp);
    
    if (argparser_keyword_present(pp, "-timeStep"))
      { o->timeStep = argparser_get_next_double(pp, 0.001, 1000.0); }
    else
      { o->timeStep = 1.0; }
    
    if (argparser_keyword_present(pp, "-groupStats"))
      { o->groupStats.tLo = (nmsim_time_t)argparser_get_next_int(pp, 0, nmsim_time_MAX); 
        o->groupStats.tHi = (nmsim_time_t)argparser_get_next_int(pp, o->groupStats.tLo, nmsim_time_MAX); 
      }
    else
      { o->groupStats.tLo = nmsim_time_MAX;  o->groupStats.tHi = 0; }
    
    argparser_get_keyword(pp, "-nSteps");
    o->nSteps = argparser_get_next_int(pp, 0, nmsim_step_count_MAX);
    
    /* Parse the external inputs spec: */
    nens_input_spec_vec_t exInput = nens_input_spec_vec_new(100);    /* Specs of external neuron inputs. */
    int32_t nex = 0; /* Counts external input specs seen. */
    while (argparser_keyword_present(pp, "-exInput"))
      { nens_elem_neuron_time_range_t kt = parse_elem_neuron_time_range(pp);
        double I_avg_raw = argparser_get_next_double(pp, -1000.0, +1000.0);
        double I_dev_raw = argparser_get_next_double(pp, 0.0, +1000.0);
        double I_avg = I_avg_raw*o->timeStep;
        double I_dev = I_dev_raw*sqrt(o->timeStep);
        nens_input_spec_t exr = (nens_input_spec_t)
          { .kt = kt, .I_avg = I_avg, .I_dev =I_dev };
        nens_input_spec_vec_expand(&exInput, nex);
        exInput.e[nex] = exr;
        nex++;
      }
    nens_input_spec_vec_trim(&exInput, nex);
    o->exInput = exInput;
    
    /* Parse the neuron trace specs: */
    nens_trace_spec_vec_t trspec = nens_trace_spec_vec_new(100);
    int32_t ntr = 0; /* Counts external input specs seen. */
    while (argparser_keyword_present(pp, "-trace"))
      { nens_elem_neuron_time_range_t kt = parse_elem_neuron_time_range(pp);
        nens_trace_spec_t trr = (nens_trace_spec_t) { .kt = kt };
        nens_trace_spec_vec_expand(&trspec, ntr);
        trspec.e[ntr] = trr;
        ntr++;
      }
    nens_trace_spec_vec_trim(&trspec, ntr);
    o->trace = trspec;
    
    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next_non_keyword(pp);

    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }

nens_elem_neuron_time_range_t parse_elem_neuron_time_range(argparser_t *pp)
  { 
    int64_t ineLo = argparser_get_next_int(pp, 0, nmsim_elem_neuron_ix_MAX);
    int64_t ineHi = argparser_get_next_int(pp, ineLo, nmsim_elem_neuron_ix_MAX);
    int64_t tLo = argparser_get_next_int(pp, 0, nmsim_time_MAX);
    int64_t tHi = argparser_get_next_int(pp, tLo, nmsim_time_MAX);
    nens_elem_neuron_time_range_t kt = (nens_elem_neuron_time_range_t)
      { .ineLo = (nmsim_elem_neuron_ix_t)ineLo,
        .ineHi = (nmsim_elem_neuron_ix_t)ineHi,
        .tLo = (nmsim_time_t)tLo,
        .tHi = (nmsim_time_t)tHi
      };
    return kt;
  }
