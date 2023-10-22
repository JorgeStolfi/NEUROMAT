#! /bin/bash 
# Last edited on 2013-12-02 13:33:31 by stolfilocal

eeg_file="$1"; shift
out_prefix="$1"; shift

PROG="nmeeg_comp_analysis"
PROGDIR=".."

${PROGDIR}/${PROG} \
  -pattern B0 data/s014_blinks_P000_eig.txt \
  -pattern B1 data/s014_blinks_P001_eig.txt \
  -pattern B2 data/s014_blinks_P002_eig.txt \
  -normalize \
  -writeComp BL0 ${out_prefix}_BL0.txt = B0 \
  -writeComp BL1 ${out_prefix}_BL1.txt = B1 \
  -writeComp BL2 ${out_prefix}_BL2.txt = B2 \
  -writeComp BL ${out_prefix}_BL.txt = B0 B1 B2 \
  < ${eeg_file} \
  > ${out_prefix}.txt

#  -norm BN = B0 B1 \
