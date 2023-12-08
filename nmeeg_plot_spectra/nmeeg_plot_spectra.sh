#! /bin/bash
# Last edited on 2023-12-05 05:10:22 by stolfi

# Usage:
#
#   nmeeg_plot_spectra.sh \
#     {SHOW} \
#     {FILE_NAME} \
#     {FMAX} \
#     {FIST_ELECTR} {LAST_ELECTR} \
#     {VSIZE} {STYLE} \
#     [ -outPrefix {OPREF} ]
#
# where
#
#   {SHOW} "SHOW" to display, "NOSHOW" not to.
#   {FILE_NAME} is the file name prefix.
#   {FMAX} is the max frequency to plot in Hz.
#   {FIRST_ELECTR} is the index of the first electrode to plot (from 0)
#   {LAST_ELECTR} is the index of the last electrode to plot (from 0), or 9999.. for last electrode.
#   {VSIZE} is the plot's total vertical size in pixels.
#   {STYLE} is a gnuplot plotting style (eg. "linespoints lt 1 pt 7") or "popsicles".
#
# The optional "-outPrefix" key specifies the output file name prefix {OPREF} (minus the ".png")
# The default {OPREF} is {FILE_NAME}.
#
# Reads {FILE_NAME}.txt, which is supposed to contain the power spectra of
# the electrodes, one frequency per line.  Writes {OPREF}.png, the
# spectral plot.

echo "$@" 1>&2

show=$1; shift
fname=$1; shift
fmax=$(echo "$1" | sed -e 's:^0*\([0-9]\):\1:'); shift;
inie=$(echo "$1" | sed -e 's:^0*\([0-9]\):\1:'); shift
fine=$(echo "$1" | sed -e 's:^0*\([0-9]\):\1:'); shift
vsize=$(echo "$1" | sed -e 's:^0*\([0-9]\):\1:'); shift; 
style="$1"; shift

opref="${fname}"
while [[ $# -gt 0 ]]; do
  if [[ "/$1" == "/-outPrefix" ]]; then
    if [[ $# -ge 2 ]]; then
      shift;
      opref="$1"; shift
    else
      echo "** missing '-outPrefix' argument" 1>&2; exit 1
    fi
  else
    echo "** spurious arguments '$*'" 1>&2; exit 1
  fi
done

if [[ ${inie} -lt 0 ]]; then
  echo "** bad first electrode index" 1>&2 ; exit 1
fi

tmp=/tmp/$$

if [[ ! ( -s ${fname}.txt ) ]]; then
  echo "** ${fname}.txt not found" 1>&2 ; exit 1
fi

ne=( `cat ${fname}.txt | egrep -m 1 -e '^ne[ ]+[=]' | sed -e 's:^.*= *0*\([0-9]\):\1:'` )
echo "found ${ne} electrodes" 1>&2 

if [[ ${fine} -ge ${ne} ]]; then
  echo "plotting all electrodes in input file" 1>&2 
  fine=$(( ${ne} - 1 ))
fi

if [[ ( ${fine} -lt ${inie} ) || ( ${fine} -ge ${ne} ) ]]; then
  echo "** bad electrode range %{inie}..${fine} should be in 0..$(( ${ne} - 1 ))" 1>&2 ; exit 1
fi

if [[ $(( ${fine} - ${inie} + 1 )) -gt 140 ]]; then
  echo "** too many electrodes" 1>&2 ; exit 1
fi

# Prepare the data file:
cat ${fname}.txt \
  | expand \
  | tr -d '\015' \
  | egrep -v -e '^[ ]*[a-zA-Z\#]' \
  > ${tmp}.dat

wc -l ${tmp}.dat
# cat ${tmp}.dat

# Get the electrode names:
ename=( `cat ${fname}.txt | grep -e '^channels =' | sed -e 's:^.*= *::'` )

# Prepare the plot command file:

color=( ee0000 aa5500 447700 008800 007755 0055ff 3333ff 5500ff )
nc=${#color[@]}

ke=$(( ${inie} ))
kc=0
sep=""
printf "plot" > ${tmp}.gpl
while [[ ${ke} -le ${fine} ]]; do
  printf "%s \\\\\n" ${sep} >> ${tmp}.gpl
  if [[ "/${style}" == "/popsicles" ]]; then 
    printf "  '%s' using 2:(rmsv(%d)) notitle with impulses lc rgb '#%s', \\\\\n" \
      "${tmp}.dat" "$(( 5 + ke ))"                    "${color[${kc}]}" >> ${tmp}.gpl
    printf "  '%s' using 2:(rmsv(%d)) title '%s' with points pt 7 lc rgb '#%s'" \
      "${tmp}.dat" "$(( 5 + ke ))"  "${ename[${ke}]}" "${color[${kc}]}" >> ${tmp}.gpl
  else
    printf "  '%s' using 2:(rmsv(%d)) title '%s' with %s lc rgb '#%s'" \
      "${tmp}.dat" "$(( 5 + ke ))" "${ename[${ke}]}" "${style}" "${color[${kc}]}" >> ${tmp}.gpl
  fi
  ke=$(( ke + 1 ))
  kc=$(( kc + 1 ))
  if [[ ${kc} -ge ${nc} ]]; then kc=0; fi
  sep=','
done
printf "\n" >> ${tmp}.gpl

# cat ${tmp}.gpl 1>&2

export GDFONTPATH=.

gnuplot <<EOF
vsize=${vsize}
set term png size 2800,(2*vsize) font "arial,16"
set output "${tmp}.png"
fmax=${fmax}
fbreak=20; # Plot f-scale breakpoint; do not change.
fstep=(fmax > 120 ? 20 : 10) # Step in second plot half - multiple of fbreak!
ymin=1.0e-3
ymax=100
rmsv(k) = sqrt(column(k))

set multiplot layout 1,2 title "RMS amplitude (uV)"
# ----------------------------------------------------------------------
# First half: 0 to 35 Hz
set origin 0.000, 0.000
set size 0.600,0.950

set xrange [-0.20:35.2]
set yrange [(0.99*ymin):(1.01*ymax)]
set logscale y

set xtics out 0,1,35
set mxtics 5
set grid xtics lt 1 lw 3 lc rgb '#dddddd', lt 1 lw 1.5 lc rgb '#dddddd'
set grid mxtics

unset y2tics
set ytics out nomirror
set mytics
set grid ytics lt 1 lw 3 lc rgb '#dddddd', lt 1 lw 1.5 lc rgb '#dddddd'
set grid mytics

set nokey
load "${tmp}.gpl"
# ----------------------------------------------------------------------
# Second half: 0 to {fmax} Hz
set origin 0.580, 0.000
set size 0.400,0.950
set key center rmargin

set xrange [-(0.01*fmax):(1.01*fmax)]
set y2range [(0.99*ymin):(1.01*ymax)]
set logscale y2

set xtics out 0,20
set mxtics 4
set grid xtics lt 1 lw 3 lc rgb '#dddddd', lt 1 lw 1.5 lc rgb '#dddddd'
set grid mxtics

set ytics out right nomirror format ""

unset ytics
set y2tics out nomirror format ''
set my2tics
set grid y2tics lt 1 lw 3 lc rgb '#dddddd', lt 1 lw 1.5 lc rgb '#dddddd'
set grid my2tics

load "${tmp}.gpl"
# ----------------------------------------------------------------------
unset multiplot
quit
EOF

plotfile="${opref}.png"

convert ${tmp}.png -resize '50%' ${plotfile}
ls -l ${plotfile}

if [[ "${show}" == "SHOW" ]]; then
  if [[ -s ${plotfile} ]]; then
    display -title '%d/%f' ${plotfile}
  else
    echo "** spectrum gnuplot failed" 1>&2; exit 1
  fi
fi

rm -fv ${tmp}.*
ls -l ${plotfile}
