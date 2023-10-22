#! /bin/bash
# Last edited on 2021-08-31 07:39:18 by stolfi

# Usage:
#
#   nmeeg_plot_channels.sh {SHOW} {FILE_NAME} {FSMP} {VMAX} {FIST_ELECTR} {LAST_ELECTR} {NSKIP} {NPLOT} {VSIZE}
#
# where
#
#   {SHOW} "SHOW" to display, "NOSHOW" not to.
#   {FILE_NAME} is the file name prefix.
#   {FSMP} is the sampling rate (samples per second) if not specified in file.
#   {VMAX} is the maximum abs value, for Y scaling purposes.
#   {FIRST_ELECTR} is the index of the first channel to plot (from 0).
#   {LAST_ELECTR} is the index of the last channel to plot (from 0) or -1 for all.
#   {NSKIP} Number of data frames to skip from beginning of file.
#   {NPLOT} Number of data frames to plot, or 0 for `until end-of-file'.
#   {VSIZE} vertical plot size in pixels.
#
# Reads {FILE_NAME}.txt, which is supposed to contain the readings of
# the channels, one time frame per line.  Writes {FILE_NAME}.png, the
# spectral plot.

echo "$@" 1>&2

show=$1; shift                # "SHOW" or "NOSHOW".
fname="$1"; shift             # EEG file name sans ".txt".
fsmp=$1; shift                # Sampling frew (0 = get from file header).
vmax=$1; shift                # Max value in vertical plot axis.
inie=$(( 10#$1 )); shift      # Index of first channel to plot (from 0).
fine=$(( 10#$1 )); shift      # INdex of last channel to plot (from 0) or -1 for last channel.
nskip=$(( 10#$1 )); shift     # Number of data frames to skip.
nplot=$(( 10#$1 )); shift     # Number of data frames to plot (-1 = to EOF).
vsize=$(( 10#$1 )); shift     # Vertical plot size in pixels.

if [[ ${inie} -lt 0 ]]; then
  echo "** bad first signal index (${inie})" 1>&2 ; exit 1
fi

if [[ ${nskip} -lt 0 ]]; then
  echo "** bad nskip (${nskip})" 1>&2 ; exit 1
fi

if [[ ! ( -e ${fname}.txt) ]]; then
  echo "** no file \"${fname}.txt\"" 1>&2 ; exit 1
fi

tmp=/tmp/$$

# Get the channel names:
chname=( `cat ${fname}.txt | egrep -m 1 -e '^channels[ ]+[=]' | sed -e 's:^.*= *::'` )

if [[ ${fine} -lt 0 ]]; then
  echo "plotting all channels in input file" 1>&2 
  fine=$(( ${#chname[@]} - 1 ))
fi
if [[ ${fine} -lt ${inie} ]]; then
  echo "** bad signal index range (${inie} .. ${fine})" 1>&2 ; exit 1
fi
if [[ ${fine} -ge ${#chname[@]} ]]; then
  echo "** invalid channel range ${inie} .. ${fine} - only ${#chname[@]} channels\n" 1>&2 ; exit 1
fi

# Number of electrodes to plot:
ne=$(( ${fine} + 1 - ${inie} ))

# Try to get the sampling frequency from the input file header: 
fsmp_h=( `cat ${fname}.txt | egrep -m 1 -e '^fsmp[ ]+[=]'` )
if [[ ${#fsmp_h[@]} -gt 0 ]]; then
  fsmp="`echo ${fsmp_h[2]} | sed -e 's:[.][0]*$::g'`";
  echo "assuming sampling frequency ${fsmp} Hz from file header" 1>&2
else
  echo "assuming sampling frequency ${fsmp} Hz from command line" 1>&2
fi
if [[ ${fsmp} -le 0 ]]; then
  echo "** bad sampling frequency (${fsmp})" 1>&2 ; exit 1
fi

# Vertical plot size before shrinking:
tmpvsize=`echo "2*${vsize}" | bc -lq`

# Prepare the data file:
cat ${fname}.txt \
  | expand \
  | tr -d '\015' \
  | egrep -v -e '^[ ]*[a-zA-Z\#]' \
  | gawk \
      -v fsmp=${fsmp} \
      -v nskip=${nskip} \
      -v nplot=${nplot} \
      -v icmax=$(( ${fine} )) \
      ' BEGIN { t=0; } 
        (t < nskip){ t++; next; } 
        ((nplot <= 0) || (t < nskip+nplot)){
          if (NF < icmax+1) { printf "**bug %d %d\n", FNR, NF > "/dev/stderr"; exit(1); }
          printf "%8d %10.6f %s\n", t, (t+0.0)/fsmp, $0; t++;
        }
      ' \
  > ${tmp}.dat

# Get the number of frames to plot: 
if [[ ${nplot} -le 0 ]]; then
  # Get the number of frames to plot from the data file: 
  echo "plotting to the end of the input file" 1>&2
  nplot="`cat ${tmp}.dat | wc -l`"
fi
if [[ ${nplot} -le 0 ]]; then
  echo "** bad nplot (${nplot})" 1>&2 ; exit 1
fi

# Choose time axis tics:
tplot=( `echo "scale=0; (${nplot}+0.0)/${fsmp}" | bc -lq` )
if [[ ${tplot} -le 15 ]]; then
  setxtics=( set xtics out 0.5\; set mxtics 5 )
elif [[ ${tplot} -le 30 ]]; then
  setxtics=( set xtics out 1.0\; set mxtics 5 )
elif [[ ${tplot} -le 60 ]]; then
  setxtics=( set xtics out 2.0\; set mxtics 4 )
elif [[ ${tplot} -le 150 ]]; then
  setxtics=( set xtics out 5.0\; set mxtics 5 )
elif [[ ${tplot} -le 300 ]]; then
  setxtics=( set xtics out 20.0\; set mxtics 4 )
elif [[ ${tplot} -le 600 ]]; then
  setxtics=( set xtics out 50.0\; set mxtics 5 )
elif [[ ${tplot} -le 1500 ]]; then
  setxtics=( set xtics out 100.0\; set mxtics 5 )
elif [[ ${tplot} -le 3000 ]]; then
  setxtics=( set xtics out 200.0\; set mxtics 4 )
elif [[ ${tplot} -le 6000 ]]; then
  setxtics=( set xtics out 500.0\; set mxtics 5 )
elif [[ ${tplot} -le 15000 ]]; then
  setxtics=( set xtics out 2000.0\; set mxtics 4 )
elif [[ ${tplot} -le 30000 ]]; then
  setxtics=( set xtics out 5000.0\; set mxtics 5 )
elif [[ ${tplot} -le 60000 ]]; then
  setxtics=( set xtics out 10000.0\; set mxtics 5 )
elif [[ ${tplot} -le 150000 ]]; then
  setxtics=( set xtics out 20000.0\; set mxtics 4 )
elif [[ ${tplot} -le 300000 ]]; then
  setxtics=( set xtics out 50000.0\; set mxtics 5 )
elif [[ ${tplot} -le 600000 ]]; then
  setxtics=( set xtics out 100000.0\; set mxtics 5 )
elif [[ ${tplot} -le 1500000 ]]; then
  setxtics=( set xtics out 200000.0\; set mxtics 4 )
else
  echo "** plot time interval too long (${tplot})" 1>&2 ; exit 1
fi

# Suppress sample dots if there are too many frames:
if [[ ${nplot} -le 600 ]]; then
  style="linespoints"
else
  style="lines"
fi

# Suppress key if there are too many channels
if [[ ${ne} -le 30 ]]; then
  setkey=( set key center rmargin )
else
  setkey=( unset key )
fi

# Line color palette:
color=( ee0000 007755 aa5500 0055ff 447700 3333ff 008800 5500ff )
nc=${#color[@]}

# Prepare the plot command file:
ke=$(( ${inie} ))
kc=0
sep=""
printf "plot" > ${tmp}.gpl
while [[ ${ke} -le ${fine} ]]; do
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
set term png size 2800,${tmpvsize} font "arial,20"
${setxtics[@]}
set grid xtics lt 1 lw 3 lc rgb '#ffddaa', lt 1 lw 1.5 lc rgb '#ffddaa'
set grid mxtics
set output "${tmp}.png"
vmax = ${vmax}
set yrange [(-vmax):(+vmax)]
${setkey[@]}
load "${tmp}.gpl"
quit
EOF

convert ${tmp}.png -resize '50%' ${fname}.png

if [[ "/${show}" == "/SHOW" ]]; then
  display -title '%d/%f' ${fname}.png
fi

rm -fv ${tmp}.*
