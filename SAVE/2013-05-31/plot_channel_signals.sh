#! /bin/bash
# Last edited on 2013-05-31 02:52:39 by stolfilocal

# Usage:
#
#   plot_channel_signals.sh {SHOW} {FILE_NAME} {RSAMP} {VMAX} {FIST_ELECTR} {LAST_ELECTR} {NSKIP} {NPLOT}
#
# where
#
#   {SHOW} "SHOW" to display, "NOSHOW" not to.
#   {FILE_NAME} is the file name prefix.
#   {RSAMP} is the sampling rate (samples per second).
#   {VMAX} is the maximum abs value, for Y range purposes.
#   {FIRST_ELECTR} is the index of the first channel to plot (from 0).
#   {LAST_ELECTR} is the index of the last channel to plot (from 0).
#   {NSKIP} Number of data frames to skip from beginning of file.
#   {NPLOT} Number of data frames to plot.
#
# Reads {FILE_NAME}.txt, which is supposed to contain the readings of
# the channels, one time frame per line.  Writes {FILE_NAME}.png, the
# spectral plot.

show=$1; shift
fname="$1"; shift
rsamp=$1; shift
vmax=$1; shift
inie=$1; shift
fine=$1; shift
nskip=$1; shift
nplot=$1; shift

if [[ 10#${inie} -lt 0 ]]; then
  echo "bad first signal index" 1>&2 ; exit 1
fi

if [[ 10#${fine} -lt 10#${inie} ]]; then
  echo "bad signal index range" 1>&2 ; exit 1
fi

if [[ 10#${nskip} -lt 0 ]]; then
  echo "bad nskip" 1>&2 ; exit 1
fi

if [[ 10#${nplot} -lt 0 ]]; then
  echo "bad nplot" 1>&2 ; exit 1
fi

if [[ 10#${nplot} -le 1000 ]]; then
  style="linespoints"
else
  style="lines"
fi

tmp=/tmp/$$

# Prepare the data file:
cat ${fname}.txt \
  | expand \
  | tr -d '\015' \
  | egrep -v -e '^[ ]*[a-zA-Z\#]' \
  | gawk \
      -v rsamp=${rsamp} \
      -v nskip=${nskip} \
      -v nplot=${nplot} \
      ' BEGIN { t=0 ;} 
        (t < nskip){ t++; next; } 
        (t < nskip+nplot){ printf "%8d %10.6f %s\n", t, (t+0.0)/rsamp, $0; t++; }
      ' \
  > ${tmp}.dat

# Get the channel names:
chname=( `cat ${fname}.txt | grep -e '^channels =' | sed -e 's:^.*= *::'` )

# Prepare the plot command file:

color=( ee0000 aa5500 447700 008800 007755 0055ff 3333ff 5500ff )
nc=${#color[@]}

ke=$(( 10#${inie} ))
kc=0
sep=""
printf "plot" > ${tmp}.gpl
while [[ ${ke} -le 10#${fine} ]]; do
  printf "%s \\\\\n" ${sep} >> ${tmp}.gpl
  printf "  '%s' using 2:%d title '%s' with %s lc rgb '#%s'" \
    "${tmp}.dat" "$(( 3 + ke ))" "${chname[${ke}]}" "${style}" "${color[${kc}]}" >> ${tmp}.gpl
  ke=$(( ke + 1 ))
  kc=$(( kc + 1 ))
  if [[ ${kc} -ge ${nc} ]]; then kc=0; fi
  sep=','
done
printf "\n" >> ${tmp}.gpl

# cat ${tmp}.gpl 1>&2

export GDFONTPATH=.

gnuplot <<EOF
set term png size 2400,1000 font "arial,20"
set xtics out 0.5
set mxtics 5
set grid xtics lt 1 lw 3 lc rgb '#ffddaa', lt 1 lw 1.5 lc rgb '#ffddaa'
set grid mxtics
set output "${tmp}.png"
vmax = ${vmax}
set yrange [(-vmax):(+vmax)]
load "${tmp}.gpl"
quit
EOF

convert ${tmp}.png -resize '50%' ${fname}.png

if [[ "${show}" == "SHOW" ]]; then
  display ${fname}.png
fi

rm -fv ${tmp}.*
