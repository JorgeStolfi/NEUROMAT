#! /usr/bin/gawk -f
# Last edited on 2013-05-28 20:56:15 by stolfilocal

# Creates a synthetic EEG data file for testing.
# The file is written to standard output.

BEGIN \
  {
    ne = 20;
    nt = 200;
    pi = 3.1415926;

    # Channel parameters:
    split("", A); # Amplitude.      
    split("", F); # Frequency.      
    split("", P); # Initial phase.  
    split("", C); # Shift.          
    for (i = 0; i < ne; i++) 
      { zi = (i + 0.5)/(ne - 0.5); # A channel parameter between 0 and 1.
        A[i] = 10 + 30.0*zi;   # Amplitude.
        F[i] = int(zi*129)/nt; # Frequency.
        P[i] = zi;             # Initial phase.
        C[i] = 50 - 100*zi;    # Shift.
        printf "channel %3d  A = %10.5f  F = %10.5f  P = %10.5f  C = %10.5f\n", i, A[i], F[i], P[i], C[i] > "/dev/stderr";
      }
    
    for (t = 0; t < nt; t++)
      { for (i = 0; i < ne; i++) 
          { val = A[i]*sin(2*pi*(F[i]*t + P[i])) + C[i];
            printf " %12.5e", val;
          }
        trig = ((t % 500) == 0 ? 200.0 : 0.00 );
        printf " %12.5e\n", trig;
      }
  }
