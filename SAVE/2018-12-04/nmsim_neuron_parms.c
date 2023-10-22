/* See {nmsim_neuron_parms.h} */
/* Last edited on 2017-07-24 21:21:11 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <nget.h>
#include <fget.h>
#include <filefmt.h>
#include <jsmath.h>

#include <nmsim_firing_func.h>
#include <nmsim_neuron_parms.h>
 
nmsim_neuron_parms_t* nmsim_neuron_parms_new
  ( double V_B,   
    double V_R,   
    double c_B,   
    double M_R,   
    double M_mu,  
    double G_R,   
    double G_mu,  
    double H_R,   
    double H_mu,  
    struct nmsim_firing_func_t *Phi
  )
  {
    nmsim_neuron_parms_t* pm = notnull(malloc(sizeof(nmsim_neuron_parms_t)), "no mem");
    (*pm) = (nmsim_neuron_parms_t)
      { .V_B = V_B, .V_R = V_R, 
        .c_B = c_B, .M_R = M_R, .M_mu = M_mu,
        .G_R = G_R, .G_mu = G_mu,
        .H_R = H_R, .H_mu = H_mu,
        .Phi = Phi
      };
    return pm;
  }

nmsim_neuron_parms_t *nmsim_neuron_parms_read(FILE *rd, double timeStep)
  {
    auto double read_param(char *name, double vmin, double vmax);
      /* Reads a line \"{name} = {value}\", including the end-of-line,
        and checks that {value} is in the range {[vmin_vmax]}. */

    auto double read_mu_from_tau_param(char *name);
      /* Reads a line \"{name} = {value}\", including the end-of-line,
       and then converts {value} from a characteristic time
       to a decay factor based on the given {timeStep}. */
    
    double read_param(char *name, double vmin, double vmax)
      { double v = nget_double(rd, name);
        if ((v < vmin) || (v > vmax))
          { fprintf
              ( stderr, "** parameter {%s} = %24.16e is out of range [ %24.16e _ %24.16e ]\n", 
                name, v, vmin, vmax
              ); 
            demand(FALSE, "aborted");
          }
        fget_eol(rd);
        return v;
      }
    
    double read_mu_from_tau_param(char *name)
      { double tau = read_param(name, 0.0, INF);
        if (tau < 0.02*timeStep)
          { return 0.0; }
        else if (tau == INF)
          { return 1.0; }
        else
          { return exp(-timeStep/tau); } 
      }

    /* Read header line: */
    filefmt_read_header(rd, nmsim_neuron_parms_FILE_TYPE, nmsim_neuron_parms_VERSION);
    
    /* Read the parameters: */
    double V_B  = read_param("V_B", -200.0, +200.0);
    double V_R  = read_param("V_R", -200.0, +200.0);
    double c_B  = read_mu_from_tau_param("tau_V");
    double M_R  = read_param("M_R", 0.0, 1.0);
    double M_mu = read_mu_from_tau_param("tau_M");
    double G_R  = read_param("G_R", 0.0, 1.0);
    double G_mu = read_mu_from_tau_param("tau_G");
    double H_R  = read_param("H_R", 0.0, 1.0);
    double H_mu = read_mu_from_tau_param("tau_H");
    
    /* Read name and parameters of the firing function: */
    char *Phi_name = nget_string(rd, "Phi");
    double V_M = fget_double(rd);
    double D_M = fget_double(rd);
    int32_t deg = fget_int32(rd);
    fget_eol(rd);
    struct nmsim_firing_func_t *Phi = nmsim_firing_func_new(Phi_name, V_M, D_M, deg);
    
    nmsim_neuron_parms_t *parms = nmsim_neuron_parms_new
      ( V_B, V_R, c_B, M_R, M_mu, G_R, G_mu, H_R, H_mu, Phi );

    /* Read footer line: */
    filefmt_read_footer(rd, nmsim_neuron_parms_FILE_TYPE);
    
    return parms;
  }
