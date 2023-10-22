#! /usr/bin/gawk -f
# Last edited on 2013-06-11 22:48:17 by stolfilocal

# Client must define {ne,nc,fsmp}

BEGIN { 
  split("", val); 
  nt = 0;
  mt = int(0.6*fsmp);
}

/^[ ]*[-+0-9]/ { 
  for (ic = 0; ic < nc; ic++)
    { val[nt,ic] = $(ic+1); }
  nt++
  next;
}

//{ 
  printf "%s:%d: ** error \n",FILENAME, FNR; 
  exit(1);
}

END {
  printf "processed %d frames, softening %d at each end\n", nt, mt > "/dev/stderr";
  for (it = 0; it < nt; it++)
    { if (it < mt) { w0 = weight(it, mt); } else { w0 = 1; }
      jt = nt - 1 - it;
      if (jt < mt) { w1 = weight(jt, mt); } else { w1 = 1; }
      for (ie = 0;  ie < ne; ie++) { printf " %20.12e", w0*w1*val[it,ie]; }
      for (ic = ne; ic < nc; ic++) { printf " %20.12e", val[it,ie]; }
      printf "\n";
    }
}

function weight(it,mt,   u)
  { u = (it + 0.5)/mt;
    u = u*u*(3 - 2*u);
    u = u*u*(3 - 2*u);
    return u;
  }
