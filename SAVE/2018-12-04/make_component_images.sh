#! /bin/bash
# Last edited on 2013-11-30 06:12:24 by stolfilocal

# Usage:
#
#   make_component_images.sh {SHOW} {FNAME}.txt ..
#
# where
#
#   {SHOW}          Either "SHOW" to display, "NOSHOW" not to.
#   {FNAME}.txt ..  One or more input single-frame EEG datasets
#
# Reads each of the listed files {FNAME}.txt.
# Assumes each file is supposed to contain header lines and 
# one EEG dataframe with values of
# EEG channels.  Writes {FNAME}.png, the image with those readings
# interpolated over the schematic head.

show=$1; shift
files=( "$@" )

set -o xtrace

btype=1 # Basis type

nshown=0
for f in ${files[@]} ; do
  cpref="${f%.*}"
  nmeeg_animate \
    ${cpref} 0 1 1 ${btype} 560 640 \
    < ${f}
  frames=( `ls ${cpref}_f000000.png | sort` )
  if [[ ${#frames[@]} -ne 1 ]]; then 
    echo "nmeeg_animate did not yield a single frame" 1>&2 ; exit 1
  fi
  mv -v "${frames[0]}" "${cpref}.png"
  if [[ ( "/${show}" == "/SHOW" ) && ( -e ${cpref}.png ) && ( ${nshown} -lt 3 ) ]]; then
    display -title '%f' ${cpref}.png
    nshown=$(( ${nshown} + 1 ))
  fi
done
