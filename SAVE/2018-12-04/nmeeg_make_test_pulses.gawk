#! /usr/bin/gawk -f
# Last edited on 2013-11-21 01:39:41 by stolfilocal

# Creates a synthetic EEG data file with various pulses for testing.
# If {type} is 0, generates only impulses at various locations.
# If {type} is 1, generates one impulse, some Gaussian humps, square pulse, chirp. 
# If {type} is 2, generates polynomial ramps of various degrees. 
# The file is written to standard output.

BEGIN \
  {
    if (nt == "") { arg_error(("must define {nt}")); }
    if (type == "") { arg_error(("must define {type}")); }

    # Assumed parameters:
    fsmp = 600;

    pi = 3.1415926;

    # Channel synthesis parameters:
    # Pulse index {ip} in each group ranges in {0..np-1} where {np} is num pulses in group.
    
    ie_imp = -1;       # Channel of first impulse, or -1.
    np_imp = 0;        # Number of impulses. 
    split("", it_imp); # Frames of impulses (index {0..np_imp}).

    ie_gauss = -1;       # Channel of first Gaussian pulse, or -1.
    np_gauss = 0;        # Number of Gaussian pulses.
    split("", it_gauss); # Central frame of pulse {ip} is {it_gauss[ip]}.
    split("", sd_gauss); # Standard deviation of pulse (in frames) is {sd_gauss[ip]}.

    ie_rect = -1; # Channel of rect pulse, or -1.
    
    ie_chirp = -1; # Channel of chirp, or -1.
    
    ie_ramp = -1; # Channel of first ramp, or -1.
    np_ramp = 0;  # Number of ramps.
    split("", dg_ramp); # Degree of ramp is {dg_ramp[ip]}.
    split("", cf_ramp); # Coeffs of ramp are {cf_ramp[ip,0..dg_ramp[ip]]}.

    ne = 0;
    if (type == 0)
      {        
        # No chirp:
        ie_chirp = -1;
        
        # Several impulses:
        ie_imp = ne;
        np_imp = 7;
        for (ip = 0; ip < np_imp; ip++)
          { fp = ip/(np_imp - 1);
            it_imp[ip] = int(2 + (nt - 5)*fp + 0.5); # Frame of impulse.
            ne++;
          }
         
        # No Gaussian pulses: 
        ie_gauss = -1;
        np_gauss = 0;

        # No rectangular pulse:
        ie_rect = -1;
        
        # No ramps:
        ie_ramp = -1;
        np_ramp = 0;
      }
    else if (type == 1)
      { 
        # One chirp:
        ie_chirp = ne; # Channel of chirp.
        it_ini_chirp = int(0.1*nt); # First frame of chirp.
        it_fin_chirp = int(0.9*nt); # Last frame of chirp.
        fr_ini_chirp = 10*fsmp/nt;  # Initial freq of chirp.
        fr_fin_chirp = 0.250*fsmp;  # Final freq of chirp.
        printf "chirp freq from %.6f to %.6f\n", fr_ini_chirp, fr_fin_chirp > "/dev/stderr";
        ne++;

        # One impulse:
        ie_imp = ne;
        np_imp = 1;
        it_imp[0] = int(0.1*nt); # Frame of impulse.
        ne++;

        # Some Gaussian pulses:
        ie_gauss = ne; #First channel of Gaussian pulses.
        np_gauss = 3; # Number of Gaussian pulses.
        sd_gauss_ini = 0.05*nt;  # SD of first Gaussian pulse.
        sd_gauss_fin = 3.00;     # SD of last Gaussian pulse.
        it_free = int(0.30*nt);  # Next "free" frame.
        for (ip = 0; ip < np_gauss; ip++)
          { fp = ip/(np_gauss-1);
            sd_gauss[ip] = exp((1-fp)*log(sd_gauss_ini) + fp*log(sd_gauss_fin));
            it_gauss[ip] = it_free + 3*sd_gauss[ip];
            it_free = it_gauss[ip] + 4*sd_gauss[ip];
            ne++;
          }

        # One rectangular pulse:
        ie_rect = ne;
        it_rect = int(0.8*nt);  # Center frame of rectangular pulse.
        ht_rect = int(0.05*nt); # Half-width of rectangular pulse.
        ne++;
        
        # No ramps:
        ie_ramp = -1;
        np_ramp = 0;
      }
    else if (type == 2)
      { 
        # No chirp:
        ie_chirp = -1;
        
        # No impulses:
        ie_imp = -1;

        # No Gaussian pulses: 
        ie_gauss = -1;
        np_gauss = 0;

        # No rectangular pulse:
        ie_rect = -1;
        
        # Several ramps:
        ie_ramp = ne;
        np_ramp = 4;
        for (ip = 0; ip < np_ramp; ip++)
          { dg = ip; # Degree of polynomial.
            split("", cf); # Coeffs of polynomial.
            # Ramp coeffs assume arg is in {[-1 _ +1]}:
            # Put {dg_ramp[ip]} roots in that interval.
            cf[0] = 1.0;
            for (k = 0; k < dg; k++)
              { ang = pi*(k + 0.5)/dg;
                rt = cos(ang);
                # Multipy polynomial by {x - rt}:
                cf[k+1] = cf[k];
                for (j = k; j > 0; j--) { cf[j] = cf[j-1] - rt*cf[j]; }
                cf[0] = - rt * cf[0];
              }
            # Save in table:
            dg_ramp[ip] = dg;
            for (j = 0; j <= dg; j++) { cf_ramp[ip,j] = cf[j]; }

            printf "channel %d ramp %d degree %d =", ne, ip, dg_ramp[ip] > "/dev/stderr";
            for (j = 0; j <= dg_ramp[ip]; j++) { printf " %+12.6f", cf_ramp[ip,j] > "/dev/stderr"; }
            printf "\n" > "/dev/stderr";
            
            ne++;
         }
      }

    # Add some trigger/marker channels:
    nc = ne + 3;
    ng = nc - ne; # Number of trigger pulses in file.

    # Write a sketchy header:
    printf "nt = %d\n", nt;
    printf "ne = %d\n", ne;
    printf "nc = %d\n", nc;
    printf "channels =";
    for (ie = 0; ie < ne; ie++) { printf " C%d", ie; }
    for (ic = ne; ic < nc; ic++) { printf " T%d", (ic - ne); }
    printf "\n";
    printf "fsmp = %.10f\n", fsmp;
    printf "orig.sample_range = %d %d\n", 0, nt-1;

    # Write the data frames:
    for (it = 0; it < nt; it++)
      { # Add the simulated electrode pulses:
        for (ie = 0; ie < ne; ie++)
          { if ((ie >= ie_imp) && (ie < ie_imp + np_imp))
              { ip = ie - ie_imp;
                if (ip >= np_imp) { printf "duh?" > "/dev/stderr"; exit(1); }
                # Single-frame impulse:
                val = (it == it_imp[ip] ? 1000.0 : 0);
              }
            else if (ie == ie_rect)
              { # Rectang pulse:
                val = ((it >= it_rect - ht_rect) && (it <= it_rect + ht_rect) ? 100.0 : 0);
              }
            else if (ie == ie_chirp)
              { # Sweep-freq sinusoid: 
                if ((it < it_ini_chirp) || (it > it_fin_chirp))
                  { val = 0; }
                else
                  { nt_chirp = it_fin_chirp - it_ini_chirp + 1;
                    u = (it - it_ini_chirp + 0.5)/nt_chirp;
                    C = log(fr_fin_chirp/fr_ini_chirp);
                    phase = (fr_ini_chirp*nt_chirp/fsmp/C)*exp(C*u);
                    f = (fr_ini_chirp*nt_chirp/fsmp)*exp(C*u);
                    v = (u < 0.5 ? u : 1-u)/0.05; # Rel. position in shoulder, from 0 end.
                    w = (v > 1 ? 1 : 0.5*(1 - cos(pi*v)));
                    val = 40*w*cos(2*pi*phase);
                    # printf "t = %10.6f phase = %10.6f freq = %10.6f val = %10.6f\n", it/fsmp, phase, f, val > "/dev/stderr";
                  }
                #Extra shift:
                val = val - 40;
              }
            else if ((ie >= ie_gauss) && (ie < ie_gauss + np_gauss))
              { # Gaussian pulse.
                ip = ie - ie_gauss;
                if (ip >= np_gauss) { printf "duh?" > "/dev/stderr"; exit(1); }
                zp = (it - it_gauss[ip])/sd_gauss[ip];
                if (zp < 0) { zp = -zp; }
                val = (zp < 8.5 ? 100*exp(-0.5*zp*zp) : 0.0);
              }
            else if ((ie >= ie_ramp) && (ie < ie_ramp + np_ramp))
              { # Polynomial ramps:
                ip = ie - ie_ramp;
                dg = dg_ramp[ip]; # Ramp degree.
                xt = 2*(it + 0.5)/nt - 1; # Frame index mapped to {[-1 _ +1]}.
                val = cf_ramp[ip,dg]; # Polynomial value.
                for (j = dg-1; j >= 0; j--) { val = xt*val + cf_ramp[ip,j]; }
                val = val*100;
              }
            val = val + 10*(ie - 0.5*(ne - 1));
            printf " %12.5e", val;
          }
        # Add the trigger pulses:
        rstage = ng*(it + 0.5)/nt + 0.5; # Stage re trigger pulses, in [0.5 _ ng+0.5]}.
        kstage = int(rstage);     # Integer part of stage, {[0..ng]}.
        fstage = rstage - kstage; # Fraction part of stage, {[0 _ 1]}
        if (kstage >= ng) { kstage = 0; }
        for (ic = ne; ic < nc; ic++)
          { jstage = (ic - ne);
            val = (10 + 4*jstage)*((kstage == jstage) && (fstage < 0.05));
            val = val - 130 - 10*(ic - ne);
            printf " %12.5e", val;
          }
        printf "\n", val;
      }
  }

function arg_error(msg) 
  { printf "** %s\n", msg > "/dev/stderr"; 
    exit(1);
  }
