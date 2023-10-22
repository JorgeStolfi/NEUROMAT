#define PROG_NAME "nmsim_mu0_diagram"
#define PROG_DESC "Generates a phase diagram of a stochastic perceptron network"
#define PROG_VERS "1.0"

#define nmsim_mu0_diagram_C_COPYRIGHT \
  "Copyright © 2016 by the State University of Campinas (UNICAMP)"

/* Last edited on 2016-05-30 21:57:13 by stolfilocal */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -Phi {PHI_NAME} {V_MID} {V_DELTA} \\\n" \
  "    -step {W_STEP} {RHO_STEP} \\\n" \
  "    -Vplot {V_MIN} {V_MAX} \\\n" \
  "    -Wplot {W_MIN} {W_MAX} \\\n" \
  "    " argparser_help_info_HELP " \\\n" \
  "    < {INFILE} \\\n" \
  "    > {OUTFILE}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program outputs an Encapsulated Postscript (EPS) file with the {(W,\rho)} phase diagram for the mean-field analysis of an all-to-all random networks of perceptron-like stochastic neurons with a specified firing function {\Phi}.\n" \
  "\n" \
  "  The neurons in the network evolve synchronously in discrete time steps.  At each time step, each neuron computes an /input potential/ {V} that is the average of its input signals received in the previous step, times a /total input weight/ {W}.  Each neuron then independently generates a boolean /output signal/ that is 1 with probability {\Phi(V)}.  When the output is 1, we say that the neuron /fired/ in that time step.\n" \
  "\n" \
  "  The firing function {\Phi} is a member of a family with name {PHI_NAME}, e.g. \"gauss\" for the Gaussian integral.  The member within that family is selected by two parameters: the /midpoint potential/ parameters {V_MID} (the membrane potential for which {\Phi} has value 1/2), and the /half-width/ {V_DELTA} (whose meaning is specific to the family}.\n" \
  "\n" \
  "  For this program, each neuron is assumed to receive inputs from all other neurons, with the same weight.  For a sufficiently large network, the state of the network just after some step {t} can be summarized by {\rho[t]}, the fraction of neurons that fired in that time step.  In the next step, all neurons will have the same input potential {V[t+1]=W\rho[t]}.  Therefore, the fraction that will fire in the next step is\n" \
  "\n" \
  "    { \rho[t+1] = \Phi(W \rho[t]) }         (1).\n" \
  "\n" \
  "  Therefore, the network is a dynamic system on the real unit interval {U=[0_1]}.  The {(W,\rho)} phase diagram shows the ultimate fate of the network when the total input weight of every neuron is {W} and the firing fraction at the end of step 0 is {\rho}.\n" \
  "\n" \
  "  If {\Phi} is sigmoidal (monotonically increasing, surjective on {U}, continuous, twice differentiable, with a single inflection point) the evolution is fairly simple.  For each {W}, the recurrence (1) has between 1 and 3 fixed points {\rho_k(W)}, all positive and less than 1.  These fixed points divide {U} into 2 to 4 sub-intervals.   For all starting states {\rho[0]} in each sub-interval, the network monotonically converges to one of the fixed points on the boundary of that interval.\n" \
  "\n" \
  "  The {(W,\rho)} phase diagram therefore has lines corresponding to the stationary states {\rho_k(W)}, and vertical lines at critical values of {W} where stationary regimes appear and disappear, and the convergence pattern changes.  These lines divide the strip {(0 _ +\oo)×[0 _ 1]} into regions of similar behavior.  The program plots that diagram, clipped to some maximum {W}.\n" \
  "\n" \
  "  If {\Phi} is not sigmoidal, the phase diagram can be fairly complex.  The program is designed to cope with firing functions that are at least monotonically non-decreasing.  For other firing functions {\Phi}, the network can have oscillating or chaotic regimes; the program is not designed to cope with them.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -PhiName {PHI_NAME} {V_MID} {V_DELTA}\n" \
  "    This mandatory argument specifies the name {PHI_NAME} of the family of the firing function {\Phi}, and the two parameters {V_MID,V_DELTA} define the function within that family.\n" \
  "\n" \
  "  -step {W_STEP} {RHO_STEP}\n" \
  "    This two mandatory argument defines the basic increments in {W} and {\rho} used to systematically explore the phase diagram, looking for stationary states.  Finer steps may be used by the program when it detects a stationary state or critical weight {W}.\n" \
  "\n" \
  "  -Vplot {V_MIN} {V_MAX}\n" \
  "    This mandatory argument specifies the range of potentials to be show in the plot of {\Phi}.\n" \
  "\n" \
  "  -Wplot {W_MIN} {W_MAX}\n" \
  "    This mandatory argument specifies the range of input weights {W} to be shown in the phase diagram.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  salamic(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2016-05-28 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  " 2016-05-28 by J. Stolfi: Created.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " nmsim_mu0_diagram_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <jsfile.h>
#include <jsmath.h>
#include <argparser.h>
#include <r2.h>
#include <pswr.h>

/* COMMAND-LINE OPTIONS */

typedef struct nsd_options_t
  { char *Phi_name;       /* Name of firing function family, e.g. "gauss". */
    double Phi_Vmid;      /* Potential where {\Phi} is 1/2. */
    double Phi_Vdelta;    /* Half-width of {\Phi} */
    double step_W;        /* Step in {W} for phase diagram scan. */
    double step_rho;      /* Step in {\rho} for phase diagram scan. */
    double Vplot_min;     /* Max {V} to show in plot of {\Phi}. */
    double Vplot_max;     /* Min {V} to show in plot of {\Phi}. */
    double Wplot_min;     /* Max {W} to show in phase diagram. */
    double Wplot_max;     /* Min {W} to show in phase diagram. */
  } nsd_options_t;

typedef struct nsd_diagram_t
  { nsd_locus_vec_t *locus;   /* Loci (lines) of stable states. */
    nsd_region_vec_t *region; /* Regions of convergence. */  
    nsd_crit_vec_t *crit;     /* Critical points. */  
  } nsd_diagram_t;

/* INTERNAL PROTOTYPES */

nsd_options_t *nsd_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and packs them as an {nsd_options_t}. */
  
int32_t main(int32_t argc,char** argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char** argv)
  {
    nsd_options_t *o = nsd_parse_options(argc, argv);
    
    nsd_diagram_t *diag = nsd_compute_diagram(o, ...);

    set term png size 2400,1000 font "courbd,24"

    double mmpt = (72.0/25.4); /* Conversion factor to turn mm into pt. */
    char *prefix = "x";
    bool_t eps = TRUE;
    char *docName = "y";
    char paperSize = NULL;
    double hSize = 6.0 * 72.0; 
    double vSize = 2.5 * 72.0; 
    PSStream *ps = pswr_new_stream("x", stdout, TRUE, "y", NULL, hSize, vSize);               }
    /* pswr_new_picture(ps, xMin, xMax, yMin, yMax);  */
    nsd_plot_Phi(ps, ..., diag);
    nsd_plot_diagram(ps, ..., diag);
    pswr_close_stream(ps);



    return 0;
  }

nsd_diagram_t nsd_compute_diagram
  ( nsd_phi_t Phi, 
    double Vmid, 
    double Vdelta, 
    double step_W, 
    double step_rho, 
    double Wmax
  )
  {
    nsd_diagram_t *diag = ndsd_new_diagram();
    int32_t nD_locus = 0;  /* Number of lines of stationary states in {diag}. */
    int32_t nD_region = 0; /* Number of regions in {diag}. */
    int32_t nD_crit = 0;   /* Number of critical points in {diag}. */
    
    /* Data for previous value of {W}: */
    double W_prev = NAN;       /* The previous value of {W}. */
    int32_t nR_prev = -1;      /* Number of stationary states, or {-1} initially. */
    double rho_prev[MAX_STAT]; /* The stationary states were {rho_prev[0..nR_prev-1]} */
    
    /* Data for the current value of {W}: */
    double W = 0.0;       /* Total input gain of neuron. */
    int32_t nR = -1;      /* Number of stationary states, or {-1} initially. */
    double rho[MAX_STAT]; /* The stationary states were {rho_prev[0..nR_prev-1]} */
    
    while (W <= Wmax) 
      { 
        nsd_compute_stable_rho_values
          ( Phi, Vmid, Vdelta,
            W, step_rho,
            &nR, rho
          );
        
        nsd_append_stable_rho_values
          ( Phi, Vmid, Vdelta,
            W_prev, nR_prev, rho_prev,
            W, nR, rho,
            diag, &nD_locus, &nD_region, &nD_crit
          );
        W = nsd_increment_parm(W, step_W, Wmax);
      }
      
    nsd_locus_vec_trim (&diag->locus,  nD_locus );
    nsd_region_vec_trim(&diag->region, nD_region);
    nsd_crit_vec_trim  (&diag->crit,   nD_crit  );
    
    return duag;
  }
  
nsd_diag_t nsd_diag_new(void)
  { 
    nsd_diag_t *diag = notnull(malloc(sizeof(nsd_diagram_t)), "no mem");
    diag->locus = nsd_locus_vec_new(100);   /* Stationary state lines. */
    diag->region = nsd_region_vec_new(100); /* Regions of convergence. */  
    diag->crit = nsd_crit_vec_new(100);     /* Critical points. */
    return diag;
  }
    

void nsd_plot_Phi(ps, ...)
  { void pswr_set_window
      ( PSStream *ps,
        double xMin, double xMax,
        double yMin, double yMax,

        double hMin, double hMax,
        double vMin, double vMax
      )
    nsd_plot_Phi_values(ps, ...);
    nsd_plot_Phi_lines(ps, ...);
    nsd_plot_Phi_frame(ps, ...);

    # Create file "${phifile}" with {(V,\Phi(V))} pairs:
    echo "Creating plot of \\Phi..." 1>&2 
    phifile="${tmp}_phi.dat"
    compute_phi_values.gawk \
        -f Phi_${PhiName}.gawk \
        -v Vref=${PhiVref} \
        -v Pref=${PhiPref} \
        -v Parm=${PhiParm} \
        -v Vmin=${Vmin} \
        -v Vmax=${Vmax} \
      > ${phifile}
  
    set origin 0.01, 0.00
    set size 0.58,1.00

    set title "Firing function Phi"

    # Command to plot the line {\rho = V/W} on the {\Phi} graph:
    printf "${sep} "'\\'"\n" >> ${phi_gplfile}
    printf "  '${W_line_file}' using 3:2 with lines title '%s'" "${Wtitle}"  >> ${phi_gplfile}
    printf " lt 1 lw 1.75 lc rgb '${color[$ic]}', "'\\'"\n" >> ${phi_gplfile}
    printf "  '${W_dots_file}' using 3:2 with points notitle"  >> ${phi_gplfile}
    printf " lt 1 lw 1.75  pt 7 ps 1.5 lc rgb '#003399'" >> ${phi_gplfile}

      # Writes to standard ouptut the pairs {(V,\Phi(V))}
      # for {V} spanning {[Vmin _ Vmax]}.

      # User must load the definition of {Phi} with 
      # -f "Phi_{NAME}.gawk"

      Parm = get_num_arg("Parm", Parm, -1000, +1000);   # Shape parameter for {\Phi}
      Vref = get_num_arg("Vref", Vref, +0.1, +99.9);    # Reference potential for {\Phi}.
      Pref = get_num_arg("Pref", Pref, 0.0000, 1.0000); # Reference probability for {\Phi}.

      Vmin = get_num_arg("Vmin", Vmin, -99.99, -0.01);
      Vmax = get_num_arg("Vmax", Vmax, Vref + 0.1, 10*Vref);

      NS = 2000;

      for (i = 0; i <= NS; i++)
        { s = (i + 0.0)/NS;
          V = Vmin + s*(Vmax - Vmin)
          rho = Phi(V,Vref,Pref,Parm);
          printf "%6.2f %7.5f\n", V, rho;
        }


    unset xrange;  set xrange  [Vmin:Vmax];
    unset x2range; set x2range [Vmin:Vmax];
    unset yrange;  set yrange  [-0.025:+1.025];
    unset y2range; set y2range [-0.025:+1.025];

    unset x2tics; unset mx2tics;
    set xtics mirror 10; set mxtics 2;
    unset format x2; set format x "%3.0f";
    set xlabel "Total input V"

    unset y2tics; unset my2tics;
    set ytics mirror 0.1; set mytics 1;
    set format y "%3.2f"; unset format y2;
    set ylabel "Firing probability"

    set grid xtics lt 1 lw 3 lc rgb '#ffddaa', lt 1 lw 1.5 lc rgb '#ffddaa'
    set grid ytics lt 1 lw 3 lc rgb '#ffddaa', lt 1 lw 1.5 lc rgb '#ffddaa'

    set key bottom right

    phi_gplfile="${tmp}_phi_plot_cmd.gpl" # Plot command for the {(V,\Phi(V))} graph.
    printf "plot "'\\'"\n" > ${phi_gplfile}
    printf "  '${phifile}' using 1:2 with lines notitle" >> ${phi_gplfile}
    printf " lt 1 lw 2.75 lc rgb '#008800'" >> ${phi_gplfile}
    
  }

void nsd_plot_diagram(ps, ...)
  { void pswr_set_window
      ( PSStream *ps,
        double xMin, double xMax,
        double yMin, double yMax,

        double hMin, double hMax,
        double vMin, double vMax
      )
    nsd_plot_diagram_regions(ps, ...);
    nsd_plot_diagram_stationary(ps, ...);
    nsd_plot_diagram_critical(ps, ...);
    nsd_plot_diagram_lines(ps, ...);
    nsd_plot_diagram_frame(ps, ...);

    Wtitle=`printf "W = %6.2f" "${Wi}"`

    # Command to plot the line {G = W/(Vref/Pref)} on the phase diagram:
    printf "${sep} "'\\'"\n" >> ${str_gplfile}
    printf "  '${W_line_file}' using 1:2 with lines title '%s'" "${Wtitle}"  >> ${str_gplfile}
    printf " lt 1 lw 1.75 lc rgb '${color[$ic]}', "'\\'"\n" >> ${str_gplfile}
    printf "  '${W_dots_file}' using 1:2 with points notitle"  >> ${str_gplfile}
    printf " lt 1 lw 1.75  pt 7 ps 1.5 lc rgb '#003399'" >> ${str_gplfile}

    set origin 0.61, 0.00
    set size 0.38,1.00

    set title "(W,rho[0]) phase diagram" 

    unset xrange;  set xrange  [Wmin_phd-0.05:Wmax_phd+0.05];
    unset x2range; set x2range [Wmin_phd-0.05:Wmax_phd+0.05];
    unset yrange;  set yrange  [-0.025:+1.025];
    unset y2range; set y2range [-0.025:+1.025];

    unset x2tics; unset mx2tics;
    set xtics mirror 25; set mxtics 5;
    unset format x2; set format x "%3.0f";
    set xlabel  "Total input weight W"

    unset ytics; unset mytics;
    set y2tics mirror 0.1; set my2tics 2;
    unset format y; set format y2 "%3.2f";
    set ylabel "rho[0]"

    set grid xtics lt 1 lw 3 lc rgb '#ffddaa', lt 1 lw 1.5 lc rgb '#ffddaa'
    set grid y2tics lt 1 lw 3 lc rgb '#ffddaa', lt 1 lw 1.5 lc rgb '#ffddaa'

    unset key

    str_gplfile="${tmp}_str_plot_cmd.gpl" # Plot command for the {(G,\rho[0])} phase diagram.
    printf "plot "'\\'"\n" > ${str_gplfile}
    printf "  '${rhofile}' using 1:2 with points notitle" >> ${str_gplfile}
    printf " lt 1 lw 2.25 pt 7 ps 0.40 lc rgb '#008800' "'\\'"\n" >> ${str_gplfile}


  }


color=( '#ff0000' '#dd4400' '#886600' '#668800' '#008800' '#0066ff' '#0022ff' '#4400ff' '#8800dd' '#dd0077' '#ff0077' )
nc=${#color[@]}   # Number of colors.


void nsd_search_roots(r0,r1,eps,tau,  rk,Vk,fk,rprev,fprev,rm)
  { 
    # Search for roots in the closed interval {[r0 _ r1]}.
    # Global parameters: {W,Vpref,Pref,Parm}.
    
    rk = r0;
    while (rk <= r1)
      { # Compute the target function at {\rho = rk}:
        Vk = rk*W;
        fk = Phi(Vk,Vref,Pref,Parm) - rk;
        if (abs(fk) <= tau)
          { # Consider {rk} a root:
            output_root(rk);
          }
        else if ((rk > r0) && (abs(fprev) > tau) && (fprev*fk < 0))
          { # We skipped over a root:
            rm = (rprev + rk)/2;
            output_root(rm);
          }
        rprev = rk;
        fprev = fk;
        # Increment {rk}:
        if (rk == r1)
          { rk += eps; }
        else if ((r1 - rk) < 1.001*eps)
          { rk = r1; }
        else if ((r1 - rk) < 2.00*eps)
          { rk = (rk + r1)/2; }
        else
          { rk += eps; }
      }
  }      

void output_root(rho,  V)
  { # Outputs a root of the equation {\rho = \Phi(\rho W)}.
    # Global parameters: {W,Vpref,Pref,Parm}.
    V = rho * W;
    # printf "(%7.5f)", rho > "/dev/stderr";
    printf "%11.5f %7.5f %6.2f\n", W, rho, V;
  }

}

void nsd_plot_special_weights()
  { 

  # Writes to standard ouptut the triples {(W,\rho,V)} 
  # for the endpoints of the {\rho = V/W} line.

  # User must load the definition of {Phi} with 
  # -f "Phi_{NAME}.gawk"
  
  Parm = get_num_arg("Parm", Parm, -1000, +1000);   # Shape parameter for {\Phi}
  Vref = get_num_arg("Vref", Vref, +0.1, +99.9);    # Reference potential for {\Phi}.
  Pref = get_num_arg("Pref", Pref, 0.0000, 1.0000); # Reference probability for {\Phi}.

  W = get_num_arg("W", W, 0.0001*Vref, 1000.0*Vref);
  
  for (i = 0; i <= 1; i++)
    { rho = i + 0.0;
      V = rho * W;
      printf "%11.5f %7.5f %6.2f\n", W, rho, V;
    }
}

#define nsd_eps_MIN (0.001)
  /* Minimum value of fundamental unit {eps} (mm).
    Must be such that coordinates that differ by this 
    much are printed differently by {nsd_write_r3_triangle}. */

#define nsd_eps_MAX (1.00)
  /* Maximum value of fundamental unit {eps} (mm). Note that good
    triangles are smaller than 1 voxel, so {eps} that is larger than
    {0.05} or so can be used only if the {step} is large too. */

#define nsd_step_MIN (0.001)
  /* Minimum value of voxel size (mm). */

#define nsd_step_MAX (200.00)
  /* Maximum value of voxel size (mm). */

nsd_options_t *nsd_parse_options(int32_t argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    nsd_options_t *o = (nsd_options_t *)malloc(sizeof(nsd_options_t)); 
    
    /* Parse keyword parameters: */
    
    /* Debug tetrahedra option: */
    o->showTetras = argparser_keyword_present(pp, "-showTetras");
    
    /* Actual voxel dimensions: */
    if (argparser_keyword_present(pp, "-step"))
      { o->step = argparser_get_next_double(pp, nsd_step_MIN, nsd_step_MAX); }
    else
      { o->step = 1.0; }

    /* Fundamental unit: */
    argparser_get_keyword(pp, "-eps");
    o->eps = (float)argparser_get_next_double(pp, nsd_eps_MIN, nsd_eps_MAX);

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);

    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }
