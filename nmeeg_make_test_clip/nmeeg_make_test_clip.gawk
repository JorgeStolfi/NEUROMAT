#! /usr/bin/gawk -f
# Last edited on 2023-11-29 15:13:43 by stolfi

# Creates a a test EEG data file by extracting selected frames and channels
# from the input real EEG file.
# The file is written to standard output.

BEGIN \
  { # The caller  must define (with "-v") the variables 
    #   {skip}    number of franes to skip.
    #   {take}    number of frames to copy.
    #   {chans}   a blank-separated list of chennel names to extract.
    #   {mkscale} (optional) a scaling factor for non-electrode channels.
    #   {mkshift} (optional) a value to be added to non-electrode channels after scaling.
    #   {mkstep}  (optional) extra value to be added to each additional non-electrode channel.
    
    # The {chans} subset must be in the same order as in the input file.
  
    if (skip == "") { arg_error(("must define {skip}")); }
    if (take == "") { arg_error(("must define {take}")); }
    if (chans == "") { arg_error(("must define {chans}")); }
    if (mkscale == "") { mkscale = 1.0; } else { mkscale = mkscale + 0.0; }
    if (mkshift == "") { mkshift = 0.0; } else { mkshift = mkshift + 0.0; }
    if (mkstep == "") { mkstep = 0.0; } else { mkstep = mkstep + 0.0; }

    nt_ot = take; 
    nc_ot = split(chans, chname_ot); # The channels are {chname_ot[1..nc_ot]}. 
    
    nc_in = -1;  # Count of channes in input file.
    ne_in = -1;  # Count of electrodes in input file.
    nt_in = -1;  # Count of frames in input file.
    split("", chname_in);  # Input channel names are {chname_in[1..nc_in]}.
    nt_read = 0; # Number of frames read. 
    nt_written = 0; # Number of frames written. 
    
    split("", val);  # The samples, indexed {[1..nt_ot,1..nc_ot]}.
    state = 0; # State 0 = expecting heading, 1 = processing frames.
 }
 
# Header processing: */

/^[a-zA-Z][_A-Za-z0-9]+ * =/ \
  { if (state != 0) { data_error("spurious header line in frames"); }
    # Keep processing.
  }

/^nt *= */ \
  { nt_in = $3 + 0; 
    if (skip + nt_ot > nt_in) { data_error(("only " nt_in " frames in input file")); }
    printf "  nt_in = %d\n", nt_in > "/dev/stderr"
    next;
  }

/^nc *= */ \
  { nc_in = $3 + 0; 
    if (nc > nc_in) { data_error(("only " nc_in " channels in input file")); }
    printf "  nc_in = %d\n", nc_in > "/dev/stderr"
    next;
  }

/^ne *= */ \
  { ne_in = $3 + 0; 
    if (ne_in > nc_in) { data_error(("inconsistent {nc} and {ne} in input file")); }
    printf "  ne_in = %d\n", ne_in > "/dev/stderr"
    next;
  }

/^channels *= */ \
  { if (NF != nc_in + 2) { data_error(("inconsistent {nc_in} = " nc_in " and channel name list")); }
    for (ic_in = 1; ic_in <= nc_in; ic_in++)
      { chname_in[ic_in] = $(ic_in+2);
        printf " %s", chname_in[ic_in] > "/dev/stderr";
      }
    printf "\n" > "/dev/stderr"
    next;
  }
  
/^fsmp *= */ \
  { fsmp = $3 + 0; 
    printf "  fsmp = %d\n", fsmp > "/dev/stderr"
    next;
  }

# Discard other headers:
/^[a-zA-Z][_A-Za-z0-9]+ * =/ \
  { next; }
  
/^ *[-+]?[0-9.]/ \
  { 
    if (state == 0)
      { # Write the header:
        printf "nt = %d\n", nt_ot;
        printf "nc = %d\n", nc_ot;
        printf "  nc_ot = %d\n", nc_ot > "/dev/stderr"
        # Whrite channel names {chname[1..nc]}and convert to {chix_ot[1..nc]}:
        split("", chix_ot);    # The 1-based indices of channels to copy are {chix_ot[1..nc]}
        ne_ot = 0; # Number of electrodes out.
        printf "channels = ";
        for (ic_ot = 1; ic_ot <= nc_ot; ic_ot++)
          { printf " %s", chname_ot[ic_ot];
            chix_ot[ic_ot] = -1;
            for (ic_in = 1; ic_in <= nc_in; ic_in++)
              { if (chname_ot[ic_ot] == chname_in[ic_in]) { chix_ot[ic_ot] = ic_in; } }
            if (chix_ot[ic_ot] < 1) { arg_error(("cannot find channel " chname_ot[ic_ot] " in input")); }
            if (chix_ot[ic_ot] <= ne_in) 
              { # Channel is an electrode:
                ne_ot++;
                if (ne_ot != ic_ot) { arg_error(("channel " chname_ot[ic_ot] " out of order")); }
              }
          }
        printf "\n";
        printf "ne = %d\n", ne_ot;
        printf "fsmp = %d\n", fsmp;
        printf "  ne_ot = %d\n", ne_ot > "/dev/stderr"
        state = 1;
      }
    nt_frames_read++;
    if (nt_frames_read < skip) { next; }
    it_ot = nt_written + 1;
    if (nt_written > nt_ot) { exit(0); }
    for (ic_ot = 1; ic_ot <= nc_ot; ic_ot++)
      { val_tc = $(chix_ot[ic_ot]);
        val[it_ot,ic_ot] = val_tc
      }
    nt_written++;
    next;
  }
  
END \
  { if (nt_written < nt_ot) { data_error("unexpected end of file"); }
    for (ic_ot = 1; ic_ot <= nc_ot; ic_ot++)
      { if (ic_ot <= ne_ot)
          { # Compute channel mean {avg} excluding outliers:
            dev = 1.0e20;
            avg = 0;
            for (pass = 1; pass <= 5; pass++)
              { # Compute channel mean {avg} excluding {avg,dev} outliers:
                sum = 0;
                num = 0;
                for (it_ot = 1; it_ot <= nt_ot; it_ot++)
                  { val_tc = val[it_ot,ic_ot];
                    if (abs(val_tc - avg) <= 4*dev) { sum += val_tc; num++; }
                  }
                if (num == 0) { prog_error("average did not converge"); }
                sum2 = 0;
                for (it_ot = 1; it_ot <= nt_ot; it_ot++)
                  { val_tc = val[it_ot,ic_ot];
                    if (abs(val_tc - avg) <= 4*dev)  { d = val_tc - avg; sum2 += d*d; }
                  }
                # Update {avg,dev}:
                avg = sum/num;
                dev = sqrt(sum2/num);
                printf "  pass %d avg = %16.7f dev = %16.7f\n", pass, avg, dev > "/dev/stderr";
              }
            # Shift channel to zero mean: 
            for (it_ot = 1; it_ot <= nt_ot; it_ot++) 
              { val[it_ot,ic_ot] = val[it_ot,ic_ot] - avg; }
          }
        else
          { # Scale marker channel values:
            for (it_ot = 1; it_ot <= nt_ot; it_ot++) 
              { val[it_ot,ic_ot] = val[it_ot,ic_ot] * mkscale + mkshift + mkstep*(ic_ot - ne_ot - 1); }
          }
      }
    
    # Write the frames:
    for (it_ot = 1; it_ot <= nt_ot; it_ot++)
      { for (ic_ot = 1; ic_ot <= nc_ot; ic_ot++)
          { printf " %16.7f", val[it_ot,ic_ot]; }
        printf "\n";
      }
  }
  
function abs(x) 
  { return (x >= 0 ? x : -x); }
  
function data_error(msg) 
  { printf "%s:%d: ** %s\n", FILENAME, FNR, msg > "/dev/stderr"; 
    printf "  [[%s]]\n", $0;
    exit(1);
  }

function arg_error(msg) 
  { printf "** %s\n", msg > "/dev/stderr"; 
    exit(1);
  }

function prog_error(msg) 
  { printf "** PROG ERROR: %s\n", msg > "/dev/stderr"; 
    exit(1);
  }
