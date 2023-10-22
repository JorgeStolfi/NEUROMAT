#! /usr/bin/gawk -f
# Last edited on 2013-11-13 21:02:50 by stolfilocal

# Creates a synthetic EEG data file for testing.
# The file is written to standard output.

BEGIN \
  {
    if (nt == "") { arg_error(("must define {nt}")); }

    # Assumed parameters:
    fsmp = 600;
    ne = 6;
    nc = 9;
    
    np = nc - ne; # Number of trigger pulses in file.

    pi = 3.1415926;

    # Channel synthesis parameters:
    fmin = 1.0; # Min frequency except for channels 0,1.
    fmax = int(nt/2) + 0.0; # Just below Nyquist's frequency. 
    split("", A); # Amplitude.      
    split("", F); # Frequency.      
    split("", P); # Initial phase.  
    split("", C); # Shift.          
    for (ie = 0; ie < ne; ie++) 
      { if (ie == 0)
          { # Channel 0 is cosine with frequency 1/2, to test end jump effect. 
            F[ie] = 0.5; P[ie] = 0.25;
          }
        else if (ie == 1)
          { # Channel 1 is sine with frequency 1/2, to test end kink effect. 
            F[ie] = 0.5; P[ie] = 0.00;
          }
        else
          { # Some sinusoid with integer frequency:
            zi = (ie - 2.0)/(ne - 3.0); # A channel parameter between 0 and 1.
            F[ie] = int(fmin*exp(zi*log(fmax/fmin)) + 0.5);     # Frequency.
            df = 2*F[ie] - nt;
            if (sqrt(df*df) < 0.00001)
              { P[ie] = 0.25; }
            eslse
              { P[ie] = zi; }
          }
        A[ie] = 200.0/sqrt(sqrt(F[ie]));
        C[ie] =  50*(ie - 0.5*(ne-1));
        printf "channel %3d  A = %10.5f  F = %10.5f  P = %10.5f  C = %10.5f\n", ie, A[ie], F[ie], P[ie], C[ie] > "/dev/stderr";
      }
    
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
      { rstage = np*(it + 0.5)/nt + 0.5; # Stage re pulses, in [0.5 _ np+0.5]}.
        kstage = int(rstage);     # Integer part of stage, {[0..np]}.
        fstage = rstage - kstage; # Fraction part of stage, {[0 _ 1]}
        if (kstage >= np) { kstage = 0; }
        for (ie = 0; ie < ne; ie++) 
          { val = A[ie]*sin(2*pi*(F[ie]*(it + 0.5)/nt + P[ie])) + C[ie];
            printf " %12.5e", val;
          }
        for (ic = ne; ic < nc; ic++)
          { if (ic == ne) 
              { val = 1000*((fstage <= 0.05) || (fstage >= 0.95)); }
            else
              { jstage = (ic - ne - 1);
                val = (50 + 20*jstage)*(kstage == jstage);
              }
            printf " %12.5e", val;
          }
        printf "\n", val;
      }
  }

function arg_error(msg) 
  { printf "** %s\n", msg > "/dev/stderr"; 
    exit(1);
  }
