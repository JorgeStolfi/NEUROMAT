#! /bin/bash
# Last edited on 2023-12-05 10:02:47 by stolfi

# Like nmeeg_plot_channels.sh, but plots two frame ranges side by side.

# Usage:
#
#   nmeeg_plot_bis_channels.sh \
#      {SHOW} {FILE_NAME} {FSMP} {VMAX} \
#      {FIST_ELECTR} {LAST_ELECTR} \
#      {NSKIP0} {NPLOT0} {HSIZE0} \
#      {NSKIP1} {NPLOT1} {HSIZE1} \
#      {VSIZE} \
#      [ -marker {MKINDEX} {MKSCALE} {MKSHIFT} ]..
#
# where
#
#   {SHOW} "SHOW" to display, "NOSHOW" not to.
#   {FILE_NAME} is the file name prefix.
#   {FSMP} is the sampling rate (samples per second) if not specified in file.
#   {VMAX} is the maximum abs value, for Y scaling purposes.
#   {FIRST_ELECTR} is the index of the first channel to plot (from 0).
#   {LAST_ELECTR} is the index of the last channel to plot (from 0) or 9999... for all.
#   {NSKIP0,NSKIP1} number of data frames to skip from beginning of file, in ranges 0 and 1.
#   {NPLOT0,NPLOT1} number of data frames to plot in ranges 0 and 1, or 0 for `until end-of-file'.
#   {HSIZE0,HSIZE1} horizontal plot size in pixels of each range plot.
#   {VSIZE} vertical plot size in pixels.
#
# Optional "-marker" parameters (repeat as needed):
#
#   {MKINDEX} index of a marker channel (from 0).
#   {MKSCALE} scale factor to apply to marker values.
#   {MKSHIFT} offset to add to marker values.
#
# Reads {FILE_NAME}.txt, which is supposed to contain the readings of
# the channels, one time frame per line.  Writes {FILE_NAME}.png, the
# spectral plot.

echo "$@" 1>&2

show=$1; shift                # "SHOW" or "NOSHOW".
fname="$1"; shift             # EEG file name sans ".txt".
fsmp=$1; shift                # Sampling frew (0 = get from file header).
vmax=$1; shift                # Max value in vertical plot axis.
inie=$(echo "$1" | sed -e 's:^0*\([0-9]\):\1:'); shift      # Index of first channel to plot (from 0).
fine=$(echo "$1" | sed -e 's:^0*\([0-9]\):\1:'); shift      # INdex of last channel to plot (from 0) or 9999... for last channel.
nskip0=$(echo "$1" | sed -e 's:^0*\([0-9]\):\1:'); shift    # Number of data frames to skip in range 0.
nplot0=$(echo "$1" | sed -e 's:^0*\([0-9]\):\1:'); shift    # Number of data frames to plot in range 0 (0 or 9999... for to EOF).
hsize0=$(echo "$1" | sed -e 's:^0*\([0-9]\):\1:'); shift    # Horizontal range 0 plot size in pixels.
nskip1=$(echo "$1" | sed -e 's:^0*\([0-9]\):\1:'); shift    # Number of data frames to skip in range 1.
nplot1=$(echo "$1" | sed -e 's:^0*\([0-9]\):\1:'); shift    # Number of data frames to plot in range 1 (0 or 9999... for to EOF).
hsize1=$(echo "$1" | sed -e 's:^0*\([0-9]\):\1:'); shift    # Horizontal range 1 plot size in pixels.
vsize=$(echo "$1" | sed -e 's:^0*\([0-9]\):\1:'); shift     # Vertical plot size in pixels.

mkops=()
opref="${fname}"
while [[ $# -gt 0 ]]; do
  if [[ "/$1" == "/-marker" ]]; then
    if [[ $# -ge 4 ]]; then
      mkops+=( "$1" "$2" "$3" "$4" );
      shift; shift; shift; shift
    else
      echo "** missing '-marker' arguments" 1>&2; exit 1
    fi
  elif [[ "/$1" == "/-outPrefix" ]]; then
    if [[ $# -ge 2 ]]; then
      shift;
      opref="$1"; shift
    else
      echo "** missing '-outPrefix' argument" 1>&2; exit 1
    fi
  else
    echo "** spurious arguments '$@'" 1>&2; exit 1
  fi
done

tmp="/tmp/$$"

tfiles=()
for kp in 0 1; do
  temp_plot_file="${tmp}_${kp}.png"
  if [[ ${kp} -eq 0 ]]; then
    nskip=${nskip0}; nplot=${nplot0}; hsize=${hsize0};
    showkey=0
  else
    nskip=${nskip1}; nplot=${nplot1}; hsize=${hsize1};
    hsize=1400; showkey=1
  fi
  nmeeg_plot_channels.sh \
      NOSHOW \
      ${fname} \
      0 ${vmax} \
      ${inie} ${fine} \
      ${nskip} ${nplot} \
      ${hsize} ${vsize} \
      -key ${showkey} \
      -outPrefix ${temp_plot_file/.png/} \
      ${mkops[@]} \
    || rm -f ${temp_plot_file}
  if [[ -s ${temp_plot_file} ]]; then 
    tfiles+=( ${temp_plot_file} )
  else
    echo "** nmeeg_plot_channels failed" 1>&2; exit 1
  fi
done

otfile="${opref}.png"   # Output ".png" file.
convert +append ${tfiles[@]} ${otfile}
if [[ -s ${otfile} ]]; then
  if [[ "/${show}" == "/SHOW" ]]; then
    display ${otfile}
  fi
else
  echo "** convert +append failed" 1>&2; exit 1
fi
rm -v ${tfiles[@]}

