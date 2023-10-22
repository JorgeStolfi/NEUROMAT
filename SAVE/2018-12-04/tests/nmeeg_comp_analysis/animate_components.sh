#! /bin/bash
# Last edited on 2013-06-10 23:11:51 by stolfilocal

# Usage:
#
#   animate_components.sh {SHOW} {PREFIX} {RSTEP}
#
# where
#
#   {SHOW} "SHOW" to display, "NOSHOW" not to.
#   {PREFIX} is the prefix of the input and output files.
#   {RSTEP} is the sub-sampling step, an integer.
#
# Reads all files named {PREFIX}_p{NNN}_esg.txt where {NNN} is a digit
# sequence. Assumes each file contains the readings of EEG channels, one
# time frame per line. Writes {PREFIX}_p{NNN}_esg.gif, the animation.

show=$1; shift
prefix="$1"; shift
rstep=$1; shift

set -o  xtrace

for f in ${prefix}_P[0-9][0-9][0-9]_esg.txt ; do \
  cpref="${f%.*}"
  nmeeg_animate \
    ${cpref} 0 0 ${rstep} 140 160 \
    < ${f}
  frames=( `ls ${cpref}_f??????.png | sort` )
  if [[ ${#frames[@]} -gt 0 ]]; then 
    convert \
        ${cpref}_m.png ${cpref}_m.png ${cpref}_m.png \
        ${frames[@]} -loop 2 -delay 5/100 \
        ${cpref}.gif ; \
    rm -f ${cpref}_f??????.png ${cpref}_m.png ${cpref}_b???.png
    if [[ ( "/${show}" == "/SHOW" ) && ( -e ${cpref}.gif ) ]]; then
      nomacs ${cpref}.gif
    fi
  fi
done
