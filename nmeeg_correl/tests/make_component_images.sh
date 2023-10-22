#! /bin/bash
# Last edited on 2021-08-29 04:36:04 by stolfi

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
ioPrefix="$1"; shift
files=( "$@" )

width=450
height=630
bdir="projects/neuromat/00-DATA/2013-09-04-ghislain-128/imgbasis/"
bpref="${bdir}/${width}_${height}_e128_"
bname="voronoi_i0_n0"

sop=( )
for f in ${files[@]} ; do
  sname="${f/${ioPrefix}/}"
  sname="${sname/.txt/}"
  sop+=( "-setName" ${sname} )
done

echo "${sop[@]}"
  
inDir="FOO"
maskDir="FOO"
basisDir="FOO"
outDir="out/FOO"

nmeeg_animate \
    -electrodes 129 \
    -basisDir ${basisDir} \
    -maskDir ${maskDir} \
    -marker 129 1.000 1.000 0.150 0.000 \
    -marker 130 1.000 0.000 0.850 1.000 \
    -inDir ${inDir} \
    -skip 0 -read 1 -step 1 \
    -outDir ${outDir} \
    ${sop[@]}

nshown=0
for f in ${files[@]} ; do
  cpref="${f%.*}"
  frames=( `ls ${cpref}_f*.png | sort` )
  if [[ ${#frames[@]} -ne 1 ]]; then 
    echo "nmeeg_animate did not yield a single frame" 1>&2 ; exit 1
  fi
  mv -v "${frames[0]}" "${cpref}.png"
  if [[ ( "/${show}" == "/SHOW" ) && ( -e ${cpref}.png ) && ( ${nshown} -lt 3 ) ]]; then
    display -title '%f' ${cpref}.png
    nshown=$(( ${nshown} + 1 ))
  fi
done
