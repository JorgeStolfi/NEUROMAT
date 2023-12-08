#! /bin/bash
# Last edited on 2023-11-04 09:02:13 by stolfi

# Reads the trace of a neuron in an elem_level simulation.
# Plots the potential the neuron as a function of time.

tmin=$1; shift;    # First time to plot.
tmax=$1; shift;    # Last time to plot.
fname="$1"; shift  # Trace file name 

# END COMMAND LINE PARSING
# ----------------------------------------------------------------------

tmp="/tmp/$$"

tfile="${tmp}.png"

# Extract the neuron number from the file name:
nnum="$(echo ${fname} | sed -e 's:[.]txt$::' -e 's:^.*_n0*\([0-9]\):\1:g')"

# Prefix for temporary file names
tmp="/tmp/$$"

export GDFONTPATH="."

gnuplot <<EOF
  set terminal png truecolor size 1600,1000 font "arial,18"
  set output "${tfile}"
  set lmargin 8; set rmargin 8; set tmargin 0.5; set bmargin 2.0
  set xrange [((${tmin})-0.5):((${tmax})+0.5)]
  set yrange [-90:+10]
  set xzeroaxis 
  set yzeroaxis
  set key bottom left
  set ytics 20
  set mytics 4
  set grid mytics lt 3
  set grid ytics lt 0
  set title "Neuron ${nnum}"
  
  # Zero if (column(k)) is 1, undef if it is 0:
  dot(k) = (column(k) == 0 ? 0/0 : 0)
  plot \
    "< egrep -e '^[ ]*[0-9]' ${fname}" using 1:2 title "V" with lines lt 1 lc rgb '#008800', \
    ""                                 using 1:8 title "J" with lines lt 1 lc rgb '#cc0000', \
    ""                                 using 1:(dot(6)) title "X" with points pt 7 ps 1.0 lc rgb '#996600'
EOF

# color=( "ff0000" "996600" "338800" "008800" "007755" "0033ff" "5500ff" "aa0066" "aa2200" "557700" "117700" "007722" "005588" "2222ff" "880088" )

if [[ -s ${tfile} ]]; then
  pfile="${fname%.*}.png"
  convert ${tfile} -resize '50%' ${pfile}
  display ${pfile}
  rm ${tfile}
fi
