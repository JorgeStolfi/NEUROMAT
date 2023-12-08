#! /bin/bash
# Last edited on 2023-12-07 14:43:04 by stolfi

# Usage:
#
#   nmeeg_plot_channels.sh \
#      {SHOW} {FILE_NAME} {FSMP} {VMAX} \
#      {FIST_ELECTR} {LAST_ELECTR} \
#      {NSKIP} {NPLOT} \
#      {HSIZE} {VSIZE} \
#      [ -key {SHOWKEY} ] \
#      [ -marker {MKINDEX} {MKSCALE} {MKSHIFT} ]..
#
# where
#
#   {SHOW} "SHOW" to display, "NOSHOW" not to.
#   {FILE_NAME} is the file name prefix.
#   {FSMP} is the sampling rate (samples per second) if not specified in file.
#   {VMAX} is the maximum abs value, for Y scaling purposes.
#   {FIRST_ELECTR} is the index of the first channel to plot (from 0).
#   {LAST_ELECTR} is the index of the last channel to plot (from 0) or 9999... for all.
#   {NSKIP} number of data frames to skip from beginning of file.
#   {NPLOT} number of data frames to plot, or 9999... for `until end-of-file'.
#   {HSIZE} horizontal plot size in pixels.
#   {VSIZE} vertical plot size in pixels.
#
# Optional "-key" parameter
#   {SHOWKEY} 1 to show the channel key (if not too many channels), 0 to omit it.
#
# Optional "-outPrefix" parameter
#   {OPREF} prefix for ".png" plot file. Default is {FILE_NAME}.
#
# Optional "-marker" parameters (repeat as needed):
#
#   {MKINDEX} index of a marker channel (from 0).
#   {MKSCALE} scale factor to apply to marker values.
#   {MKSHIFT} offset to add to marker values.
#
# Reads {FILE_NAME}.txt, which is supposed to contain the readings of
# the channels, one time frame per line.  Writes {FILE_NAME}.png, the
# spectral plot.

echo "$@" 1>&2

show=$1; shift                # "SHOW" or "NOSHOW".
fname="$1"; shift             # EEG file name sans ".txt".
fsmp=$1; shift                # Sampling frew (0 = get from file header).
vmax=$1; shift                # Max value in vertical plot axis.
inie=$(echo "$1" | sed -e 's:^0*\([0-9]\):\1:'); shift      # Index of first channel to plot (from 0).
fine=$(echo "$1" | sed -e 's:^0*\([0-9]\):\1:'); shift      # INdex of last channel to plot (from 0) or 9999... for last channel.
nskip=$(echo "$1" | sed -e 's:^0*\([0-9]\):\1:'); shift     # Number of data frames to skip.
nplot=$(echo "$1" | sed -e 's:^0*\([0-9]\):\1:'); shift     # Number of data frames to plot (0 or 9999... for to EOF).
hsize=$(echo "$1" | sed -e 's:^0*\([0-9]\):\1:'); shift     # Horizontal plot size in pixels.
vsize=$(echo "$1" | sed -e 's:^0*\([0-9]\):\1:'); shift     # Vertical plot size in pixels.

mkindex=()
mkscale=()
mkshift=()
showkey=1
opref="${fname}"
while [[ $# -gt 0 ]]; do
  if [[ "/$1" == "/-marker" ]]; then
    if [[ $# -ge 4 ]]; then
      shift
      mkindex+=( "$1" ); shift
      mkscale+=( "$1" ); shift
      mkshift+=( "$1" ); shift
    else
      echo "** missing '-marker' arguments" 1>&2; exit 1
    fi
  elif [[ "/$1" == "/-key" ]]; then
    if [[ $# -ge 2 ]]; then
      shift
      showkey="$1"; shift
    else
      echo "** missing '-key' argument" 1>&2; exit 1
    fi
  elif [[ "/$1" == "/-outPrefix" ]]; then
    if [[ $# -ge 2 ]]; then
      shift
      opref="$1"; shift
    else
      echo "** missing '-key' argument" 1>&2; exit 1
    fi
  else
    echo "** spurious arguments '$@'" 1>&2; exit 1
  fi
done

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
nc=${#chname[@]}

if [[ ${fine} -ge ${nc} ]]; then
  echo "plotting all channels in input file" 1>&2 
  fine=$(( ${nc} - 1 ))
fi
if [[ ${fine} -lt ${inie} ]]; then
  echo "** bad signal index range (${inie} .. ${fine})" 1>&2 ; exit 1
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

# Plot size before shrinking:
tmphsize=`echo "2*${hsize}" | bc -lq`
tmpvsize=`echo "2*${vsize}" | bc -lq`

efile="${fname}.txt"

# Prepare the data file:
cat ${efile} \
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
        ((nplot==0)||(t < nskip+nplot)){
          if (NF < icmax+1) { printf "**bug %d %d\n", FNR, NF > "/dev/stderr"; exit(1); }
          printf "%8d %10.6f %s\n", t, (t+0.0)/fsmp, $0; t++;
        }
      ' \
  > ${tmp}.dat
ndata="`cat ${tmp}.dat | wc -l`" # Number of frames in data file: 
echo "skipped ${nskip} frames kept/found ${ndata} frames" 1>&2

# Get the number of frames to plot: 
if [[ ${nplot} -gt ${ndata} ]]; then nplot=${ndata}; fi
if [[ ${nplot} -eq 0 ]]; then nplot=${ndata}; fi
if [[ ${nplot} -le 0 ]]; then
  echo "** bad nplot (${nplot})" 1>&2 ; exit 1
fi

# Compute time range to plot in msec, adjusted for {hsize}:
tplot=( `echo "scale=0; t=1400000*(${nplot}+${fsmp}-1)/${fsmp}/${hsize}; t+0" | bc -lq` )
# Choose time axis tics:
if [[ ${tplot} -le 1500 ]]; then
  setxtics=( set xtics out 0.05\; set mxtics 5 )
elif [[ ${tplot} -le 3000 ]]; then
  setxtics=( set xtics out 0.1\; set mxtics 5 )
elif [[ ${tplot} -le 6000 ]]; then
  setxtics=( set xtics out 0.2\; set mxtics 4 )
elif [[ ${tplot} -le 15000 ]]; then
  setxtics=( set xtics out 0.5\; set mxtics 5 )
elif [[ ${tplot} -le 30000 ]]; then
  setxtics=( set xtics out 1.0\; set mxtics 5 )
elif [[ ${tplot} -le 60000 ]]; then
  setxtics=( set xtics out 2.0\; set mxtics 4 )
elif [[ ${tplot} -le 150000 ]]; then
  setxtics=( set xtics out 5.0\; set mxtics 5 )
elif [[ ${tplot} -le 300000 ]]; then
  setxtics=( set xtics out 20.0\; set mxtics 4 )
elif [[ ${tplot} -le 600000 ]]; then
  setxtics=( set xtics out 50.0\; set mxtics 5 )
elif [[ ${tplot} -le 1500000 ]]; then
  setxtics=( set xtics out 100.0\; set mxtics 5 )
elif [[ ${tplot} -le 3000000 ]]; then
  setxtics=( set xtics out 200.0\; set mxtics 4 )
elif [[ ${tplot} -le 6000000 ]]; then
  setxtics=( set xtics out 500.0\; set mxtics 5 )
elif [[ ${tplot} -le 15000000 ]]; then
  setxtics=( set xtics out 2000.0\; set mxtics 4 )
elif [[ ${tplot} -le 30000000 ]]; then
  setxtics=( set xtics out 5000.0\; set mxtics 5 )
elif [[ ${tplot} -le 60000000 ]]; then
  setxtics=( set xtics out 10000.0\; set mxtics 5 )
elif [[ ${tplot} -le 150000 ]]; then
  setxtics=( set xtics out 20000.0\; set mxtics 4 )
elif [[ ${tplot} -le 300000000 ]]; then
  setxtics=( set xtics out 50000.0\; set mxtics 5 )
elif [[ ${tplot} -le 600000000 ]]; then
  setxtics=( set xtics out 100000.0\; set mxtics 5 )
elif [[ ${tplot} -le 1500000000 ]]; then
  setxtics=( set xtics out 200000.0\; set mxtics 4 )
else
  echo "** plot time interval too long (${tplot})" 1>&2 ; exit 1
fi
echo "tplot = ${tplot} '${setxtics[*]}'" 1>&2 

# Suppress sample dots if there are too many frames:
if [[ ${nplot} -le 600 ]]; then
  style="linespoints"
else
  style="lines"
fi

# Suppress key if there are too many channels

if [[ ${showkey} -eq 0 ]]; then
  setkey=( unset key )
elif [[ ${ne} -le 30 ]]; then
  setkey=( set key center rmargin )
else
  setkey=( unset key )
fi

# Line color palette:
color=( ee0000 007755 aa5500 0055ff 447700 3333ff 008800 5500ff )
nc=${#color[@]}

# Prepare the plot command file:
printf "plot" > ${tmp}.gpl
kc=0
sep=""
# Plot electrodes:
ke=$(( ${inie} ))
while [[ ${ke} -le ${fine} ]]; do
  printf "%s \\\\\n" ${sep} >> ${tmp}.gpl
  printf "  '%s' using 2:%d title '%s' with %s lc rgb '#%s'" \
    "${tmp}.dat" "$(( 3 + ke ))" "${chname[${ke}]}" "${style}" "${color[${kc}]}" >> ${tmp}.gpl
  ke=$(( ke + 1 ))
  kc=$(( kc + 1 ))
  if [[ ${kc} -ge ${nc} ]]; then kc=0; fi
  sep=','
done
# Plot markers, if any:
im=0
while [[ ${im} -lt ${#mkindex[@]} ]]; do
  km=${mkindex[${im}]}
  msc="${mkscale[${im}]}"
  msh="${mkshift[${im}]}"
  printf "%s \\\\\n" ${sep} >> ${tmp}.gpl
  printf "  '%s' using 2:(%s*column(%d)+%s) title '%s' with lines lc rgb '#%s'" \
    "${tmp}.dat" "${msc}" "$(( 3 + km ))" "${msh}" "${chname[${km}]}" "${color[${kc}]}" >> ${tmp}.gpl
  im=$(( im + 1 ))
  kc=$(( kc + 1 ))
  if [[ ${kc} -ge ${nc} ]]; then kc=0; fi
  sep=','
done
printf "\n" >> ${tmp}.gpl

# cat ${tmp}.gpl 1>&2

export GDFONTPATH=.

gnuplot <<EOF
set term png size ${tmphsize},${tmpvsize} font "arial,20"
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

pfile="${opref}.png"
convert ${tmp}.png -resize '50%' ${pfile}

if [[ "/${show}" == "/SHOW" ]]; then
  if [[ -s ${pfile} ]]; then
    display -title '%d/%f' ${pfile}
  else
    echo "** ${pfile} not created" 1>&2; exit 1
  fi
fi

rm -fv ${tmp}.*
