#! /bin/bash
# Last edited on 2023-11-29 11:25:43 by stolfi

inprefix="$1"; shift
runs="$1"; shift
outprefix="$1"; shift
show="$1"; shift
animate="$1"; shift

PROGDIR=".."
PROG="nmeeg_cleanup"

${PROGDIR}/${PROG} \
    -inPrefix "${inprefix}" \
    -runs ${runs//,/ } \
    -blinkTerms 2 \
    -blinkThresh 200.0 \
    -blinkRad 20 \
    -changeThresh 1.0 \
    -numCorr 3 \
    -corrThresh 0.75 \
    -outPrefix "${outprefix}"
    
echo "=== computing spectra and plots ==="
for ff in ${outprefix}*.txt ; do \
  prefix=${ff%.txt} ; 
  nmeeg_spectrum ${prefix} 0 0 < ${prefix}.txt ; \
  vmax=180
  nmeeg_plot_channels.sh ${show} ${prefix} 0 ${vmax}  0 9999  0 0  1400 500 ; \
  nmeeg_plot_spectra.sh ${show} ${prefix}_pwr 100 0 130 480 'linespoints pt 7'
  if [[ "/${animate}" == "/ANIMATE" ]]; then \
    echo "animating the principal component patterns" 1>&2 ; \
    # ??? animate_components.sh ${show} ${out_prefix} ${N_USE} ${RSTEP} ; \
  fi
done ; \

  
