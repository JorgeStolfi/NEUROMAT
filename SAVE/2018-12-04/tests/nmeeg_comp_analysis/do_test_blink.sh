#! /bin/bash 
# Last edited on 2013-12-04 11:06:39 by stolfilocal

# Blink pulse pattern only (assumed to be eigenvalue P000)

subjid="$1"; shift
runid="$1"; shift
eeg_file="$1"; shift
out_prefix="$1"; shift

PROG="nmeeg_comp_analysis"
PROGDIR=".."

${PROGDIR}/${PROG} \
  -pattern BL0 data/s${subjid}_blinks_P000_eig.txt \
  -normalize \
  -writeComp BL0 ${out_prefix}_BL0.txt = BL0 \
  -delete BL0 \
  < ${eeg_file} \
  > ${out_prefix}.txt
