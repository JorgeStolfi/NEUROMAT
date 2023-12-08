#! /bin/bash
# Last edited on 2023-12-05 10:52:00 by stolfi

# Usage:
#
#   nmeeg_plot_spectra.sh {SHOW} {CHNAME} {IN_NAME} {OT_NAME}
#
# where
#
#   {SHOW} "SHOW" to display, "NOSHOW" not to.
#   {CHNAME} is the electrode name, e.g. "C001" or "Oz".
#   {IN_NAME} is the input file with the histogram data for that electrode.
#   {OT_NAME} is the output file with the histogram plot.
#
# Reads {IN_NAME}.txt, which is supposed to contain the histogram
# of the potential of electrode {CHNAME}.  Writes {OT_NAME}.png, the
# plot of the same.

echo "$@" 1>&2

show="$1"; shift
chname="$1"; shift
inname="$1"; shift
otname="$1"; shift

tmp=/tmp/$$

if [[ ! ( -s ${inname}.txt ) ]]; then
  echo "** ${inname}.txt not found" 1>&2 ; exit 1
fi

tdfile="${tmp}_hist.txt"
tpfile="${tmp}_hist.png"

# Prepare the data file:
cat ${inname}.txt \
  | gawk \
      ' BEGIN { vhi_prev = -inf; vstep = -1; nh = 0 }
        /^ *[-+0-9.]/ { 
          if (NF != 3) { bug("bad NF") }
          vlo = $1; vhi = $2; ct = $3;
          vstep = vhi-vlo;
          if (vhi_prev == -inf) { print vlo - 0.0001*vstep, 0; }
          print(vlo, ct);
          print(vhi, ct);
          nh++;
          vhi_prev = vhi;
          next;
        }
        //{ bug("bad format"); }
        END{ 
          if (nh == 0) { bug("no bins"); }
          print vhi_prev + 0.0001*vstep, 0;
        }
      ' \
  > ${tdfile}

export GDFONTPATH=.

gnuplot <<EOF
set term png size 2800,1000 noenhanced font "arial,16"
set output "${tpfile}"
set yrange [-0.001:]
set nokey
set title "${inname}"
plot "${tdfile}" using 1:2 with lines lw 2 lc rgb '#008800'
quit
EOF

if [[ -s ${tpfile} ]]; then
  convert ${tpfile} -resize '50%' ${otname}.png
  if [[ "${show}" == "SHOW" ]]; then
    display -title '%d/%f' ${otname}.png
  else
    ls -l ${otname}.png 1>&2
  fi
else
  echo "** ${tpfile} not created" 1>&2 ; exit 1
fi
rm -fv ${tmp}.*
