#! /bin/bash 
# Last edited on 2023-11-29 11:25:31 by stolfi

# Blink pulse pattern only (assumed to be eigenvalue P000)

subjid="$1"; shift
runid="$1"; shift
eeg_file="$1"; shift
out_prefix="$1"; shift
show="$1"; shift
animate="$1"; shift

PROG="nmeeg_comp_analysis"
PROGDIR=".."

#   -delete BL0 \

echo "=== decomposring signals ===" 1>&2
${PROGDIR}/${PROG} \
  -pattern BL0 data/s${subjid}_blinks_P000_eig.txt \
  -normalize \
  -writeComp BL0 ${out_prefix}_BL0.txt = BL0 \
  < ${eeg_file} \
  > ${out_prefix}.txt

echo "=== computing spectra and plots ==="
for ff in ${out_prefix}{_*,}.txt ; do \
  prefix=${ff%.txt} ; 
  nmeeg_spectrum ${prefix} 0 0 < ${prefix}.txt ; \
  if [[ "/${prefix}" == "/${out_prefix}" ]]; then
    vmax=50
  else
    vmax=180
  fi
  nmeeg_plot_channels.sh ${show} ${prefix} 0 ${vmax}  0 9999  0 0  1400 500 ; \
  nmeeg_plot_spectra.sh ${show} ${prefix}_pwr 100 0 127 480 'linespoints pt 7'
  if [[ "/${animate}" == "/ANIMATE" ]]; then \
    echo "animating the principal component patterns" 1>&2 ; \
    # ??? animate_components.sh ${show} ${out_prefix} ${N_USE} ${RSTEP} ; \
  fi
done ; \

