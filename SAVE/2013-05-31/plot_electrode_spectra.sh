#! /bin/bash
# Last edited on 2013-05-31 04:51:09 by stolfilocal

# Usage:
#
#   plot_electrode_spectra.sh {SHOW}  {FILE_NAME} {FIST_ELECTR} {LAST_ELECTR}
#
# where
#
#   {SHOW} "SHOW" to display, "NOSHOW" not to.
#   {FILE_NAME} is the file name prefix.
#   {FIRST_ELECTR} is the index of the first electrode to plot (from 0)
#   {LAST_ELECTR} is the index of the last electrode to plot (from 0)
#
# Reads {FILE_NAME}.txt, which is supposed to contain the power spectra of
# the electrodes, one frequency per line.  Writes {FILE_NAME}.png, the
# spectral plot.

show=$1; shift
fname=$1; shift
inie=$1; shift
fine=$1; shift

if [[ 10#${inie} -lt 0 ]]; then
  echo "bad first electrode index" 1>&2 ; exit 1
fi

if [[ 10#${fine} -lt 10#${inie} ]]; then
  echo "bad electrode range" 1>&2 ; exit 1
fi

tmp=/tmp/$$

# Prepare the data file:
cat ${fname}.txt \
  | expand \
  | tr -d '\015' \
  | egrep -v -e '^[ ]*[a-zA-Z\#]' \
  > ${tmp}.dat

# Get the electrode names:
ename=( `cat ${fname}.txt | grep -e '^electrodes =' | sed -e 's:^.*= *::'` )

# Prepare the plot command file:

color=( ee0000 aa5500 447700 008800 007755 0055ff 3333ff 5500ff )
nc=${#color[@]}

ke=$(( 10#${inie} ))
kc=0
sep=""
printf "plot" > ${tmp}.gpl
while [[ ${ke} -le 10#${fine} ]]; do
  printf "%s \\\\\n" ${sep} >> ${tmp}.gpl
  printf "  '%s' using 2:%d title '%s' with impulses lc rgb '#%s', \\\\\n" \
    "${tmp}.dat" "$(( 5 + ke ))" "${ename[${ke}]}" "${color[${kc}]}" >> ${tmp}.gpl
  printf "  '%s' using 2:%d title '%s' with points pt 7 lc rgb '#%s'" \
    "${tmp}.dat" "$(( 5 + ke ))" "${ename[${ke}]}" "${color[${kc}]}" >> ${tmp}.gpl
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
set xtics out 10
set grid xtics lt 1 lw 3 lc rgb '#ffddaa', lt 1 lw 1.5 lc rgb '#ffddaa'
set output "${tmp}.png"
set yrange [0.000001:]
set logscale y
load "${tmp}.gpl"
quit
EOF

convert ${tmp}.png -resize '50%' ${fname}.png
ls -l ${fname}.png

if [[ "${show}" == "SHOW" ]]; then
  display ${fname}.png
fi

rm -fv ${tmp}.*
ls -l ${fname}.png
