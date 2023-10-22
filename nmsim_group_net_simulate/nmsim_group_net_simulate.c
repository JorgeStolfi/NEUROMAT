#define PROG_NAME "nmsim_group_net_simulate"
#define PROG_DESC "group-level simulation of a GL neuron network"
#define PROG_VERS "1.0"

/* Last edited on 2020-12-06 19:58:06 by jstolfi */ 

#define PROG_COPYRIGHT \
  "Copyright © 2019  State University of Campinas (UNICAMP)"
  
#define PROG_AUTH \
  "Created by J. Stolfi on 2019-05-28"
  
#define PROG_HIST
  
#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -netFile {NET_NAME} \\\n" \
  "    [ -timeStep {TIME_STEP} ] \\\n" \
  "    -nSteps {N_STEPS} \\\n" \
  "    [ -exInput {KLO} {KHI} {TLO} {THI} {I_AVG} {I_DEV} ].. \\\n" \
  "    [ -trace {KLO} {KHI} {TLO} {THI} ].. \\\n" \
  "    -outPrefix {OUT_PREFIX} \\\n" \
  "    " argparser_help_info_HELP " \\\n" \
  "    < {INFILE}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program reads a group-level description of" \
  " a Galves-Loecherbach neuron net and simulates its" \
  " behavior for a specified number of steps, with specified" \
  " external inputs.  Traces of selected neuron groups are saved" \
  " to disk files.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -netFile {NET_NAME} \n" \
  "    This mandatory command line argument specifies the" \
  " file name {NET_NAME} containing the description of the" \
  " network.\n" \
  "\n" \
  "  -outPrefix {OUT_PREFIX} \n" \
  "    This mandatory argument specifies the common" \
  " prefix {OUT_PREFIX} for all output file names.\n" \
  "\n" \
  "\n" \
  "  -timeStep {TIME_STEP} \n" \
  "    This optional argument specifies the simulation" \
  " time step length {TIME_STEP} in milliseconds.  If omitted," \
  " assumes \"-timeStep 1\" (1 millisecond).\n" \
  "\n" \
  "  -nSteps {N_STEPS} \n" \
  "    This mandatory argument specifies the" \
  " number {N_STEPS} of time steps to simulate.  The" \
  " discrete times will range from 0 to {N_STEPS}.  (However, the" \
  " simulation parameters {X,I,J} will not be defined for the last time {N_STEPS}.)\n" \
  "\n" \
  "  -exInput {KLO} {KHI} {TLO} {THI} {I_AVG} {I_DEV}\n" \
  "    This argument specifies the next-step external input voltage" \
  " increment {I[k][t]} for the neuron groups with indices {k} in" \
  " the range {KLO .. KHI}, for the discrete" \
  " times {t} in the range {TLO .. THI}. It can be" \
  " repeated as many times as desired.  The input value will" \
  " be a random variable with Gaussian" \
  " distribution, mean {I_AVG}, and deviation {I_DEV}.  All" \
  " those variables will be generated independently of each" \
  " other and of any other inputs or parameters.\n" \
  "\n" \
  "    !!! FIX THIS TOO !!!\n" \
  "\n" \
  "    If the same" \
  " pair {k,t} is specified in two or more of these" \
  " specs, the corresponding values are added together.  For any" \
  " neuron groups and times not specified via these arguments, the" \
  " external input {I[k][t]} will be zero.\n" \
  "\n" \
  "  -trace {KLO} {KHI} {TLO} {THI}\n" \
  "    This optional argument specifies that the program should write" \
  " out the state of each neuron group with index {k} in {KLO..KHI} for" \
  " each discrete time {t} in {TLO..THI}.  This argument can specified" \
  " more than once, to trace different neuron groups over different time" \
  " intervals; but the neuron group index ranges had better be disjoint.   The data" \
  " for each neuron group is written to a separate file" \
  " called \"{OUT_PREFIX}_g{III}.txt\" where {III} is the" \
  " neuron group index {k} formatted as 10 decimal digits, zero  padded.\n" \
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
#include <affirm.h>

#include <nmsim_basic.h>
#include <nmsim_class_net.h>
#include <nmsim_group_net.h>
#include <nmsim_group_net.h>
#include <nmsim_group_net_trace.h>

#include <nmsim_group_net_sim.h>

typedef struct ngns_neuron_group_time_range_t
  { nmsim_group_neuron_ix_t ilo; /* Index of first neuron group in set. */
    nmsim_group_neuron_ix_t ihi; /* Index of last neuron group in set. */
    nmsim_time_t tlo;                  /* First time in range. */
    nmsim_time_t thi;                  /* Last time in range. */
  } ngns_neuron_group_time_range_t;
  /* Specifies a range of neuron group indices {ilo..ihi} 
    and a range of times {tlo..thi}. */

typedef struct ngns_input_spec_t
  { ngns_neuron_group_time_range_t kt;  /* Range of neuron groups and times. */
    double I_avg;                       /* Mean input value. */
    double I_dev;                       /* Deviation of input values. */
  } ngns_input_spec_t;
  /* Description of external inputs for a range of neuron groups
    during a range or times.  Each input {I[k][t]} will be 
    a random variable independently drawn from a Gaussian 
    distribution with mean {I_avg} and deviation {I_dev}. */
    
vec_typedef(ngns_input_spec_vec_t,ngns_input_spec_vec,ngns_input_spec_t);

typedef struct ngns_trace_spec_t
  { 
    ngns_neuron_group_time_range_t kt;    /* Range of neuron groups and times. */
  } ngns_trace_spec_t;
  /* Specification of a range of neuron groups whose evolution
    should be monitored during a specified time range. */
    
vec_typedef(ngns_trace_spec_vec_t,ngns_trace_spec_vec,ngns_trace_spec_t);
    
typedef struct ngns_options_t
  { char *netFile;                /* Name of file with the network's description. */
    double timeStep;                  /* Step length to use in simulation (ms). */
    nmsim_step_count_t nSteps;        /* Number of steps to simulate. */
    ngns_input_spec_vec_t exInput;    /* Specs of external neuron group inputs. */
    ngns_trace_spec_vec_t trace;      /* Neuron Groups and times to trace. */
    char *outPrefix;                  /* Prefix for output file names. */
  } ngns_options_t;
  /* Arguments from command line. */
 
nmsim_group_net_t *ngns_read_network(char *fname, double  timeStep);
  /* Reads a descriptin of a group-level GL neuron group net 
    from file {rd}. */

nmsim_group_net_trace_t *ngns_make_trace
  ( nmsim_group_net_t *gnet, 
    nmsim_time_t tlo,
    nmsim_time_t thi,
    ngns_trace_spec_vec_t trspec
  );
  /* Creates a trace data structure for the neuron group net {gnet},
    that will store the full states of neuron groups and times 
    specified in {trspec}, clipped to the time range {tlo..thi}. */
    
typedef void ngns_trace_spec_visit_proc_t 
  ( nmsim_group_neuron_ix_t ilo_r, 
    nmsim_group_neuron_ix_t ihi_r, 
    nmsim_time_t tlo_r, 
    nmsim_time_t thi_r
  );
  /* Type of a procedure that processes a range of neuron group indices
    {ilo_r..ihi_r} and a range of times {tlo_r..thi_r},
    called by {ngns_scan_trace_specs} below.  The argument
    ranges will be non-empty and contained in the global 
    index and time intervals. */
    
void ngns_scan_trace_specs
  ( ngns_trace_spec_vec_t trspec,
    nmsim_group_neuron_count_t nng,
    nmsim_time_t tlo,
    nmsim_time_t thi,
    bool_t verbose,
    ngns_trace_spec_visit_proc_t *visit
  );
  /* Scans the trace specifications {trspec}.  For each entry
    whose neuron group index interval intersects the global index interval {0..nng-1} and 
    whose time interval intersects the global time interval {tlo..thi}, calls {visit}
    with the intersections of those intervals.
    
    If {verbose} is true, prints warnings to {stderr}
    when the {trspec} entry intervals are not contained in the
    global intervals. */

void ngns_set_external_inputs
  ( nmsim_group_neuron_count_t nng, 
    double Iavg[], 
    double Idev[], 
    nmsim_time_t t,
    ngns_input_spec_vec_t exInput
  );
  /* Defines the group external inputs {Iavg[k],Idev[k]} for the time
    step from {t} to {t+1}, for all neuron groups {k} in {0..nng-1},
    according to the specs in {exInput}. Namely, for entry {exr} of
    {exInput} such that {t} is in the time range {exr.tlo .. exr.thi},
    and for every neuron group index {k} in {exr.ilo .. exr.ihi}, sets
    {Iavg[k],Idev[k]} to the mean and deviation specified in
    {exInput}. If the same pair {k,t} is specified multiple times in
    {exInput}, the corresponding averages abd variances are added together. For any
    neuron groups and times not specified in {exInput}, the values
    {Iavg[k],Idev[k]} are set to zero. */

void ngns_simulate
  ( nmsim_group_net_t *gnet, 
    ngns_input_spec_vec_t exInput,
    nmsim_time_t tlo,
    nmsim_time_t thi,
    nmsim_group_net_trace_t *gtrace
  );
  /* Simulates the evolution of the network {gnet} from time {t=tlo}
    to time {t=thi} -- that is, for {thi-tlo} time steps.  Stores
    into {gtrace} the full states of some neuron groups in some time ranges,
    as specified in {gtrace}. 
    
    Note that the states for {t=thi} are only partially defined. */
  
nmsim_group_net_t *ngns_read_network(char *fname, double timeStep);
  /* Reads a network description from {fname}.  The {timeStep}
    is used to convert characteristc decay times to decay factors. */

void ngns_write_trace(char *prefix, nmsim_group_net_trace_t *gtrace);
  /* Writes the group-level traces in {gtrace} with names "{prefix}_group_n{NN}.txt"
    where {NN} is the neuron group index. */

ngns_options_t *ngns_parse_options(int argc, char **argv);
  /* Parses the command line arguments. */

ngns_neuron_group_time_range_t parse_neuron_group_time_range(argparser_t *pp);
  /* Parses four arguments {ILO} {IHI} {TLO} {THI} 
    as a range of neuron group indices {ILO..IHI} and a range
    of times {TLO..THI}. */

int main(int argc, char **argv);

/* IMPLEMENTATIONS: */

int main(int argc, char **argv)
  { 
    ngns_options_t *o = ngns_parse_options(argc, argv);
    nmsim_group_net_t *gnet = ngns_read_network(o->netFile, o->timeStep);
    nmsim_time_t tlo = 0;                 /* Time at start of first step of simulation. */
    nmsim_time_t thi = tlo + o->nSteps; /* Time at end of last step of simulation. */
    nmsim_group_net_trace_t *gtrace = ngns_make_trace(gnet, tlo, thi, o->trace);
    ngns_simulate(gnet, o->exInput, tlo, thi, gtrace);
    ngns_write_trace(o->outPrefix, gtrace);
    nmsim_group_net_trace_free(gtrace);
    return 0;
  }
    
void ngns_simulate
  ( nmsim_group_net_t *gnet, 
    ngns_input_spec_vec_t exInput,
    nmsim_time_t tlo,
    nmsim_time_t thi,
    nmsim_group_net_trace_t *gtrace
  )
  {
    nmsim_class_net_t *cnet = gnet->cnet; /* Description of neuron group and synapse classes. */
    
    nmsim_group_neuron_count_t nng = gnet->nng; /* Number of neuron groups. */
    
    /* Allocate work arrays {state,Iavg,Idev}: */
    nmsim_group_neuron_state_t *state = notnull(malloc(nng*sizeof(nmsim_group_neuron_state_t)), "no mem");
    double *Iavg = notnull(malloc(nng*sizeof(double)), "no mem");
    double *Idev = notnull(malloc(nng*sizeof(double)), "no mem");
    
    /* Initialize the neuron group ages, potentials, and modulators: */
    for (nmsim_group_neuron_ix_t ing = 0; ing < nng; ing++) 
      { state[ing] = nmsim_group_neuron_state_make(); /* Empty state.*/
        nmsim_class_neuron_ix_t inc = gnet->ngrp.e[ing].inc; /* Neuron class index. */
        nmsim_group_neuron_state_throw(&(state[ing]), cnet->nclass.e[inc], 2);
      }
    
    /* Simulate: */
    for (nmsim_time_t t = tlo; t < thi; t++)
      { /* Define the external inputs for this time step: */
        ngns_set_external_inputs(nng, Iavg, Idev, t, exInput);
        /* Apply the evolution equations: */
        nmsim_group_net_sim_step(gnet, t, Iavg, Idev, state, gtrace);
      }
    /* Save the non-step data for the last time {t=nSteps}: */
    if (gtrace != NULL) { nmsim_group_net_trace_set_state_vars(gtrace, thi, nng, state); }  
    
    /* Reclaim storage: */
    for (nmsim_group_neuron_ix_t ing = 0; ing < nng; ing++) 
      { nmsim_group_neuron_state_free(&(state[ing])); }
    free(state);
    free(Iavg);
    free(Idev);
  }
    
void ngns_scan_trace_specs
  ( ngns_trace_spec_vec_t trspec,
    nmsim_group_neuron_count_t nng,
    nmsim_time_t tlo,
    nmsim_time_t thi,
    bool_t verbose,
    ngns_trace_spec_visit_proc_t *visit
  )
  { for (nmsim_group_neuron_count_t k = 0; k < trspec.ne;k++)
      { ngns_trace_spec_t *spk = &(trspec.e[k]); /* Next neuron group range. */
        /* Get and clip the neuron group index range: */
        nmsim_group_neuron_ix_t ilo_clp, ihi_clp;
        char *iname = (verbose ? "neuron group index" : NULL);
        nmsim_int32_range_clip(spk->kt.ilo, spk->kt.ihi, 0, nng-1, &ilo_clp, &ihi_clp, iname); 
        /* Get and clip the time range: */
        nmsim_time_t tlo_clp, thi_clp;
        char *tname = (verbose ? "time" : NULL);
        nmsim_int64_range_clip(spk->kt.tlo, spk->kt.thi, tlo, thi, &tlo_clp, &thi_clp, tname); 
        /* If not empty, visit the ranges: */
        if ((tlo_clp < thi_clp) && (ilo_clp <= ihi_clp))
          { visit(ilo_clp, ihi_clp, tlo_clp, thi_clp); }
      }
  }

nmsim_group_net_trace_t *ngns_make_trace
  ( nmsim_group_net_t *gnet, 
    nmsim_time_t tlo, 
    nmsim_time_t thi,
    ngns_trace_spec_vec_t trspec
  )
  {
    nmsim_group_neuron_count_t nng = gnet->nng; /* Num of neuron groups in network. */
        
    /* Count neuron groups to be traced: */
    
    nmsim_group_neuron_count_t nng_tr = 0; /* Number of neuron groups to trace. */
        
    auto void rcount
      ( nmsim_group_neuron_ix_t ilo_r, 
        nmsim_group_neuron_ix_t ihi_r, 
        nmsim_time_t tlo_r, 
        nmsim_time_t thi_r
      );
      /* Assumes that the neuron group index range {ihi_r..ilo_r} has been clipped to {0..nng-1},
        and that the time range {tlo_r..thi_r} has been cliped to the simulation time range
        {tlo..thi}.  Increments {nng_tr} by the number of neuron groups in that range,
        that is, {ihi_r - ilo_r + 1}. */
    
    ngns_scan_trace_specs(trspec, nng, tlo, thi, TRUE, rcount);
    
    /* Create trace structure: */
    
    nmsim_group_net_trace_t *gtrace = nmsim_group_net_trace_new(tlo, thi, nng_tr);  /* To be allocated later. */
    nmsim_group_neuron_count_t tne = 0; /* Counts entries added to {gtrace}. */
        
    auto void rstore
      ( nmsim_group_neuron_ix_t ilo_r, 
        nmsim_group_neuron_ix_t ihi_r, 
        nmsim_time_t tlo_r, 
        nmsim_time_t thi_r
      );
      /* Assumes that the neuron group index range {ihi_r..ilo_r} has been clipped to {0..nng-1},
        and that the time range {tlo_r..thi_r} has been cliped to the simulation time range
        {tlo..thi}.  Adds to {gtrace} one entry for each neuron group index in {ilo_r..ih_r}
        with time range {tlo_r..thi_r}. */
    
    ngns_scan_trace_specs(trspec, nng, tlo, thi, FALSE, rstore);
    assert(tne == nng_tr);
    
    return gtrace;
    
    /* INTERNAL IMPLEMENTATIONS */
        
    void rcount
      ( nmsim_group_neuron_ix_t ilo_r, 
        nmsim_group_neuron_ix_t ihi_r, 
        nmsim_time_t tlo_r, 
        nmsim_time_t thi_r
      )
      { assert((0 <= ilo_r) && (ilo_r <= ihi_r) && (ihi_r < nng));
        assert((tlo <= tlo_r) && (tlo_r <= thi_r) && (thi_r <= thi));
        nng_tr += (ihi_r - ilo_r +1); 
      }
      
    void rstore
      ( nmsim_group_neuron_ix_t ilo_r, 
        nmsim_group_neuron_ix_t ihi_r, 
        nmsim_time_t tlo_r, 
        nmsim_time_t thi_r
      )
      { assert((0 <= ilo_r) && (ilo_r <= ihi_r) && (ihi_r < nng));
        assert((tlo <= tlo_r) && (tlo_r <= thi_r) && (thi_r <= thi));
        for (nmsim_group_neuron_ix_t ing = ilo_r; ing <= ihi_r; ing++)
          { gtrace->trng[tne] = nmsim_group_neuron_trace_new(ing, tlo_r, thi_r);
            tne++;
          }
      }
  }
        
void ngns_set_external_inputs
  ( nmsim_group_neuron_count_t nng, 
    double Iavg[], 
    double Idev[], 
    nmsim_time_t t,
    ngns_input_spec_vec_t exInput
  )
  {
    /* Clear the inputs: */
    for (nmsim_group_neuron_ix_t ing = 0; ing < nng; ing++) 
      { Iavg[ing] = 0.0; Idev[ing] = 0.0; } 
      
    /* Get the external inputs for this time step: */
    for (int32_t r = 0; r < exInput.ne; r ++)
      { /* Grab the input spec: */
        ngns_input_spec_t *exr = &(exInput.e[r]);
        if ((t >= exr->kt.tlo) && (t <= exr->kt.thi))
          { nmsim_group_neuron_ix_t ilo = (nmsim_group_neuron_ix_t)imax(0, exr->kt.ilo);
            nmsim_group_neuron_ix_t ihi = (nmsim_group_neuron_ix_t)imin(nng-1, exr->kt.ilo);
            
            for (nmsim_group_neuron_ix_t ing = ilo; ing <= ihi; ing++) 
              { Iavg[ing] += exr->I_avg;
                Idev[ing] = hypot(Idev[ing], exr->I_dev);
              }
          }
      }
  }

nmsim_group_net_t *ngns_read_network(char *fname, double timeStep)
  { FILE *rd = open_read(fname, TRUE);
    nmsim_group_net_t *gnet = nmsim_group_net_read(rd, timeStep);
    fclose(rd);
    return gnet;
  }

void ngns_write_trace(char *prefix, nmsim_group_net_trace_t *gtrace)
  { char *epref = NULL;
    asprintf(&epref, "%s_group", prefix);
    nmsim_group_net_trace_write(epref, gtrace);
    free(epref);
  }
    
vec_typeimpl(ngns_input_spec_vec_t,ngns_input_spec_vec,ngns_input_spec_t);
    
vec_typeimpl(ngns_trace_spec_vec_t,ngns_trace_spec_vec,ngns_trace_spec_t);
         
ngns_options_t *ngns_parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    /* Allocate the command line argument record: */
    ngns_options_t *o = notnull(malloc(sizeof(ngns_options_t)), "no mem");

    /* Parse keyword parameters: */

    argparser_get_keyword(pp, "-netFile");
    o->netFile = argparser_get_next_non_keyword(pp);
    
    if (argparser_keyword_present(pp, "-timeStep"))
      { o->timeStep = argparser_get_next_double(pp, 0.001, 1000.0); }
    else
      { o->timeStep = 1.0; }
    
    argparser_get_keyword(pp, "-nSteps");
    o->nSteps = argparser_get_next_int(pp, 0, nmsim_step_count_MAX);
    
    /* Parse the external inputs spec: */
    ngns_input_spec_vec_t exInput = ngns_input_spec_vec_new(100);    /* Specs of external neuron group inputs. */
    int32_t nex = 0; /* Counts external input specs seen. */
    while (argparser_keyword_present(pp, "-exInput"))
      { ngns_neuron_group_time_range_t kt = parse_neuron_group_time_range(pp);
        double I_avg = argparser_get_next_double(pp, -1000.0, +1000.0);
        double I_dev = argparser_get_next_double(pp, 0.0, +1000.0);
        ngns_input_spec_t exr = (ngns_input_spec_t)
          { .kt = kt, .I_avg = I_avg, .I_dev =I_dev };
        ngns_input_spec_vec_expand(&exInput, nex);
        exInput.e[nex] = exr;
        nex++;
      }
    ngns_input_spec_vec_trim(&exInput, nex);
    o->exInput = exInput;
    
    /* Parse the neuron group trace specs: */
    ngns_trace_spec_vec_t trspec = ngns_trace_spec_vec_new(100);
    int32_t ntr = 0; /* Counts external input specs seen. */
    while (argparser_keyword_present(pp, "-trace"))
      { ngns_neuron_group_time_range_t kt = parse_neuron_group_time_range(pp);
        ngns_trace_spec_t trr = (ngns_trace_spec_t) { .kt = kt };
        ngns_trace_spec_vec_expand(&trspec, ntr);
        trspec.e[ntr] = trr;
        ntr++;
      }
    ngns_trace_spec_vec_trim(&trspec, ntr);
    o->trace = trspec;
    
    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next_non_keyword(pp);

    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }

ngns_neuron_group_time_range_t parse_neuron_group_time_range(argparser_t *pp)
  { 
    int64_t ilo = argparser_get_next_int(pp, 0, nmsim_group_neuron_ix_MAX);
    int64_t ihi = argparser_get_next_int(pp, ilo, nmsim_group_neuron_ix_MAX);
    int64_t tlo = argparser_get_next_int(pp, 0, nmsim_time_MAX);
    int64_t thi = argparser_get_next_int(pp, tlo, nmsim_time_MAX);
    ngns_neuron_group_time_range_t kt = (ngns_neuron_group_time_range_t)
      { .ilo = (nmsim_group_neuron_ix_t)ilo,
        .ihi = (nmsim_group_neuron_ix_t)ihi,
        .tlo = (nmsim_time_t)tlo,
        .thi = (nmsim_time_t)thi
      };
    return kt;
  }
