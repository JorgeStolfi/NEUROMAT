#! /bin/bash
# Last edited on 2023-11-04 09:01:38 by stolfi

# Reads the trace of a neuron in an elem_level simulation.
# Plots the potential the neuron as a function of time.

tmin=$1; shift;    # First time to plot.
tmax=$1; shift;    # Last time to plot.
fname="$1"; shift  # Trace file name 

# END COMMAND LINE PARSING
# ----------------------------------------------------------------------

# Prefix for temporary file names
tmp="/tmp/$$"

tfile="${tmp}_J.png"
pfile="${fname%.*}_J.png"
  
# Extract the neuron number from the file name:
nnum="$(echo ${fname} | sed -e 's:[.]txt$::' -e 's:^.*_n0*\([0-9]\):\1:g')"

export GDFONTPATH="."

gnuplot <<EOF
  hpx = 1600
  vpx = 300
  tmin = ${tmin}
  tmax = ${tmax}
  tfile = "${tfile}"
  nnum = "${nnum}"
  vars = "total input {J}"
  load "plot_common.gpl"
  set yrange [-10:+15]
  set tmargin 0.5; set bmargin 0.5
  set ytics 5
  set mytics 5
  
  plot \
    "< egrep -e '^[ ]*[0-9]' ${fname}" using 1:8 title "J" with lines lt 1 lc rgb '#cc0000'
EOF

if [[ -s ${tfile} ]]; then
  convert ${tfile} -resize '50%' ${pfile}
  display ${pfile}
  rm ${tfile}
fi
