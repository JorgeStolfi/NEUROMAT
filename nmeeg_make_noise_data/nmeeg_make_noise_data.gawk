#! /usr/bin/gawk -f
# Last edited on 2021-08-30 18:20:24 by stolfi

# Creates a synthetic EEG data file with various pulses for testing.
# The file is written to standard output.

BEGIN \
  {
    pi = 3.1415926;
    fsmp = 500
    ne = 50;

    # Channel synthesis parameters:
    # Pulse index {ip} in each group ranges in {0..np-1} where {np} is num pulses in group.
    
    # Add some trigger/marker channels:
    nc = ne + 3;
    ng = nc - ne; # Number of trigger pulses in file.

    # Write a sketchy header:
    printf "nt = %d\n", nt;
    printf "ne = %d\n", ne;
    printf "nc = %d\n", nc;
    printf "channels =";
    for (ie = 0; ie < ne; ie++) { printf " C%03d", ie; }
    for (ic = ne; ic < nc; ic++) { printf " T%02d", (ic - ne); }
    printf "\n";
    printf "fsmp = %.10f\n", fsmp;
    printf "orig.sample_range = %d %d\n", 0, nt-1;
    
    amp = 15.0 # Amplitude (uV).

    # Write the data frames:
    for (it = 0; it < nt; it++)
      { for (ie = 0; ie < ne; ie++) 
          { val = amp*(rand()+rand()+rand()+rand()-2);
            printf " %12.5e", val;
          }
        # Add the trigger pulses:
        kstage = int(ng*(it + 0.5)/nt + 0.5); # Stage in {[0..ng]}.
        for (ic = ne; ic < nc; ic++)
          { ig = ic - ne;
            val = (ig == kstage ? 2*amp : 0)
            printf " %12.5e", val;
          }
        printf "\n";
      }
  }

function arg_error(msg) 
  { printf "** %s\n", msg > "/dev/stderr"; 
    exit(1);
  }
