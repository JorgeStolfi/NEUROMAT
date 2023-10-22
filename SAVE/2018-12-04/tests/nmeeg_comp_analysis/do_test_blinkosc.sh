#! /bin/bash 
# Last edited on 2013-12-04 11:06:21 by stolfilocal

# Blink pulse pattern plus 10 Hz oscillation (assumed to be eigenvalues P000 and P001)

subjid="$1"; shift
runid="$1"; shift
eeg_file="$1"; shift
out_prefix="$1"; shift

PROG="nmeeg_comp_analysis"
PROGDIR=".."

${PROGDIR}/${PROG} \
  -pattern BL0 data/s${subjid}_blinks_P000_eig.txt \
  -pattern H10 data/s${subjid}_blinks_P001_eig.txt \
  -normalize \
  -writeComp BL0 ${out_prefix}_BL0.txt = BL0 \
  -writeComp H10 ${out_prefix}_H10.txt = H10 \
  -writeComp BLH ${out_prefix}_BLH.txt = BL0 H10 \
  -delete BL0 H10 \
  < ${eeg_file} \
  > ${out_prefix}.txt

#  -norm BN = B0 B1 \
#  \
#   -pattern B2 data/s013_blinks_P002_eig.txt \
#   -writeComp BL2 ${out_prefix}_BL2.txt = B2 \
#   -writeComp BL ${out_prefix}_BL.txt = B0 B1 B2
