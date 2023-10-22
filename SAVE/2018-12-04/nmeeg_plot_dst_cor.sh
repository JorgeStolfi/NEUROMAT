#! /bin/bash
# Last edited on 2013-11-29 05:07:28 by stolfilocal

# Usage:
#
#   nmeeg_plot_dst_cor.sh {SHOW} {DATA_DIR} {OUT_DIR} {FILE_NAME}
#
# where
#
#   {SHOW} "SHOW" to display, "NOSHOW" not to.
#   {DATA_DIR} is the file name prefix.
#   {OUT_DIR} is the file name prefix.
#   {FILE_NAME} is the file name prefix.
#
# Reads "{DATA_DIR}/{FILE_NAME}.txt", which is supposed to contain 
# distance and correlation data for electrode pairs.
#
# Writes "{OUT_DIR}/{FILE_NAME}.png", the corresponding plot. The input
# file must have one line per electrode pair, with fields "{I} {J}
# {D2IJ} {CVIJ} {CVII} {CVJJ}" where {D2IJ} is the square of the
# distance between the 3D positions of electrodes {I} and {J}, {CVIJ} is
# the covariance of those elecrodes, and {CVII} and {CVJJ} are their
# variances.

show=$1; shift
datadir="$1"; shift
outdir="$1"; shift
fname="$1"; shift

datafile="${datadir}/${fname}.txt"
plotfile="${outdir}/${fname}.png"

if [[ ! ( -e ${datafile}) ]]; then
  echo "** no file \"${datafile}\"" 1>&2 ; exit 1
fi

tmp=/tmp/$$

export GDFONTPATH=.

gnuplot <<EOF
set term png size 2800,1000 font "arial,20"
set xrange [-0.01:2.01]
set xtics out 0.5
set mxtics 5
set grid xtics lt 1 lw 3 lc rgb '#ffddaa', lt 1 lw 1.5 lc rgb '#ffddaa'
set grid mxtics
set xzeroaxis lt 1 lw 3 lc rgb '#ffddaa'
set output "${tmp}.png"
vmax = 1.1
set yrange [(-vmax):(+vmax)]
unset key

dist(k) = sqrt(column(k))
correl(i,j,k) = column(i)/sqrt(column(j)*column(k))

plot "${datafile}" using (dist(3)):(correl(4,5,6)) with points pt 7 lc rgb '#0022ff'
quit
EOF

convert ${tmp}.png -resize '50%' ${plotfile}

if [[ "/${show}" == "/SHOW" ]]; then
  display ${plotfile}
fi

rm -fv ${tmp}.*
