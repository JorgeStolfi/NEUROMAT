#! /bin/bash
# Last edited on 2023-11-04 09:20:01 by stolfi

# Plots a scattergram of two selected channels from an EEG recording file.

show=$1; shift             # "SHOW" or "NOSHOW".
inFile="$1"; shift         # Input EEG file.
step="$(echo "$1" | sed -e 's:^0*\([0-9]\):\1:')"; shift # Plot only one evry this many frames.
xname="$1"; shift          # Name of channel to plot on X axis.
xmax="$1"; shift           # Max X axis value, for plot.
yname="$1"; shift          # Name of channel to plot on Y axis.
ymax="$1"; shift           # Max Y axis value, for plot.

pltSize=800 # Size of whole plot in pixels.

tmp=/tmp/$$

fname="${inFile/.txt/}"
pltFile="${fname}_${xname}_${yname}_scatter.png"

# Get the channel names:
chnames=( `cat ${inFile} | egrep -m 1 -e '^channels[ ]+[=]' | sed -e 's:^.*= *::'` )
echo "channels = ${chnames[*]}" 1>&2

# Vertical plot size before shrinking:
tmpSize=`echo "2*${pltSize}" | bc -lq`

# Prepare the data file:
cat ${inFile} \
  | expand \
  | gawk \
      -v step=${step} \
      -v xname=${xname} \
      -v yname=${yname} \
      -v chnames="${chnames[*]}" \
      ' BEGIN {
          nc = split(chnames, chn);
          icx = -1; icy = -1; # X and Y channel indices, from 1.
          for (ic = 1; ic <= nc; ic++) {
            chni = chn[ic];
            if (xname == chni) { icx = ic; }
            if (yname == chni) { icy = ic; }
          }
          if (icx < 1) { arg_error(("bad channel name \"" xname "\"")); }
          if (icy < 1) { arg_error(("bad channel name \"" yname "\"")); }
          it = 0 # Frame index, from 0.
        } 
        /^ *([#a-zA-Z]|$)/ { next; }
        /^ *[-+0-9.]/ {
          if ((it % step) == 0) {
            print $(icx), $(icy);
          }
          it++;
          next;
        }
        // { data_error("invalid format"); }
        
        function arg_error(msg) { printf "** %s\n", msg > "/dev/stderr"; exit(1); }
        function data_error(msg) { printf "%s%d** %s\n  [[%s]]\n", FILENAME, FNR, msg, $0 > "/dev/stderr"; exit(1); }
      ' \
  > ${tmp}.dat

wc -l ${tmp}.dat

export GDFONTPATH="${HOME}/tt-fonts"

gnuplot << EOF
set term png size ${tmpSize},${tmpSize} font "arial,20"

set output "${tmp}.png"

xmax = ${xmax}
set xrange [(-xmax):(+xmax)]
set xlabel "${xname}"
set xtics 50
set mxtics 5
set grid xtics lt 1 lw 3 lc rgb '#ffddaa', lt 1 lw 1.5 lc rgb '#ffddaa'
set grid mxtics

ymax = ${ymax}
set yrange [(-ymax):(+ymax)]
set ylabel "${yname}"
set ytics 50
set mytics 5
set grid ytics lt 1 lw 3 lc rgb '#ffddaa', lt 1 lw 1.5 lc rgb '#ffddaa'
set grid mytics

set nokey

plot "${tmp}.dat" using 1:2 with points pt 7 ps 1.0 lc rgb '#ff0000'
quit
EOF

if [[ -s ${tmp}.png ]]; then 
  convert ${tmp}.png -resize '50%' ${pltFile}
  ls -l ${pltFile}
  if [[ "/${show}" == "/SHOW" ]]; then
    display -title '%d/%f' ${pltFile}
  fi
  rm -fv ${tmp}.*
else
  echo "plot file not generated" 1>&2 ; exit 1
fi
