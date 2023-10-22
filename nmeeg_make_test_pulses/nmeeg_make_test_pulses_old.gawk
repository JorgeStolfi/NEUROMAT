#! /usr/bin/gawk -f


# Creates a synthetic EEG signal channel file with various pulses for testing.
# The file is written to standard output.

BEGIN \
  {
    if (nt == "") { arg_error(("must define {nt}")); }
    if (fsmp == "") { arg_error(("must define {fsmp}")); }

    # If {1} is 1, generates one impulse, some Gaussian humps, square pulse, chirp. 

    pi = 3.1415926;

    # Channel synthesis parameters:
    # Pulse index {ip} in each group ranges in {0..np-1} where {np} is num pulses in group.
    
    nc = 1; # Only one channel in file.
    
    nt_space = int(5.0*fsmp); # Pause before and after each pulse (frames).
    it_next = 2*nt_space; # Index of next available frame.

    # One impulse:
    np_imp = 1;        # Number of impulses. 
    split("", it_imp); # Frames of impulses (index {0..np_imp-1}).
    it_imp[0] = it_next + nt_space; # Frame of impulse.
    it_next = it_imp[0] + nt_space;
    printf "impulse at frame %d\n", it_imp[0] > "/dev/stderr";

    # Some Gaussian pulses:
    np_gauss = 5;        # Number of Gaussian pulses.
    split("", it_gauss); # Central frame of pulse {ip} is {it_gauss[ip]} (index {0..np_gauss-1}).
    split("", sd_gauss); # Standard deviation of pulse (in frames) is {sd_gauss[ip]}.
    sd_gauss_ini = 0.05*fsmp;  # SD of first Gaussian pulse (frames).
    sd_gauss_fin = 2.00*fsmp;  # SD of last Gaussian pulse (frames).
    for (ip = 0; ip < np_gauss; ip++)
      { fp = ip/(np_gauss-1);
        sd_gauss[ip] = exp((1-fp)*log(sd_gauss_ini) + fp*log(sd_gauss_fin));
        it_gauss[ip] = it_next + nt_space + int(6*sd_gauss[ip]);
        it_next = it_gauss[ip] + int(6*sd_gauss[ip]) + nt_space;
        printf "gaussian pulse at frame %d  deviation = %.6f\n", it_gauss[ip], sd_gauss[ip] > "/dev/stderr";
      }

    # One rectangular pulse:
    hw_rect = int(4*fsmp); # Half-width on either side of center.
    it_rect = it_next + nt_space + 2*hw_rec;  # Center frame of rectangular pulse.
    it_next = it_rect + 2*hw_rec + nt_space;
    printf "rect pulse on frames %d ± %d\n", it_rect, hw_rect > "/dev/stderr";

    # One chirp (sweep-freq sinusoid):
    it_ini_chirp = it_next + nt_space;             # First frame of chirp.
    it_fin_chirp = int(it_ini_chirp + 120*fsmp) - 1; # Last frame of chirp.
    fr_ini_chirp = 0.2;  # Initial freq of chirp (Hz).
    fr_fin_chirp = 0.250*fsmp;  # Final freq of chirp (Hz).
    printf "chirp in frames %d .. %d  freq from %.6f to %.6f\n", it_ini_chirp, it_fin_chirp, fr_ini_chirp, fr_fin_chirp > "/dev/stderr";
    it_next = it_fin_chirp + 1 + nt_space;
    
    it_next += nt_space; # One extra space at end.

    if (it_next > nt) 
      { printf "** BUG it_next = %d nt = %d\n", it_nex, nt > "/dev/stderr";
        exit(1);
      }

    # Write a sketchy header:
    printf "nt = %d\n", nt;
    printf "nc = %d\n", nc;
    printf "ne = %d\n", ne;
    printf "channels = SYN\n", ne;
    printf "fsmp = %.10f\n", fsmp;
    printf "sample_range = %d %d\n", 0, nt-1;

    # Write the data frames:
    for (it = 0; it < nt; it++)
      { vtot = 0;
        # Add the simulated pulses:
        # Compute the pulse's value {val} ref to its own zero:

        # Single-frame impulses:
        for (ip = 0; ip < np_imp; ip++)
          { ip = 0;
            val = (it == it_imp[ip] ? 500.0 : 0);
            vtot += val;
          }

        # Rectang pulse:
        val = ((it >= it_rect - hw_rect) && (it <= it_rect + hw_rect) ? 100.0 : 0);
        vtot += val;

        # Chirp:
        if ((it >= it_ini_chirp) && (it <= it_fin_chirp))
          { nt_chirp = it_fin_chirp - it_ini_chirp + 1;      # Frame count of chirp.
            u = (it - it_ini_chirp + 0.5)/nt_chirp;          # Rel position inside chirp, 0 to 1.
            C = log(fr_fin_chirp/fr_ini_chirp);              # Log of ratios of final to initial freq.
            phase = (fr_ini_chirp*nt_chirp/fsmp/C)*exp(C*u); # Phase of sinusoid.
            freq = fr_ini_chirp*exp(C*u);                    # Freq in Hz; derivative of phase rel to time.
            v = (u < 0.5 ? u : 1-u)/0.05;                    # Rel. position in shoulder, from 0 end.
            w = (v > 1 ? 1 : 0.5*(1 - cos(pi*v)));
            val = 40*w*cos(2*pi*phase);
            # printf "t = %10.6f phase = %10.6f freq = %10.6f val = %10.6f\n", it/fsmp, phase, freq, val > "/dev/stderr";
            vtot += val;
          }

        # Gaussian pulses:
        for (ip = 0; ip < np_gauss; ip++)
          { zp = (it - it_gauss[ip])/sd_gauss[ip];
            if (zp < 0) { zp = -zp; }
            val = (zp < 8.5 ? 100*exp(-0.5*zp*zp) : 0.0);
            vtot += val;
          }

        printf "%+9.3f\n", vtot
     }
  }

function arg_error(msg) 
  { printf "** %s\n", msg > "/dev/stderr"; 
    exit(1);
  }
