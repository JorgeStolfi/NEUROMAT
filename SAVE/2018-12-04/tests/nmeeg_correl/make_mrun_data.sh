#! /bin/bash
# Last edited on 2013-06-13 15:33:17 by stolfilocal

# Arguments:
subject="$1"; shift          # Subject ID (3 digits).
filter=$(( 10#$1 )); shift   # Use filtered (1) or unfiltered (0) runs?
rini=$(( 10#$1 )); shift     # First run.
rfin=$(( 10#$1 )); shift     # Last run.
otag="$1"; shift;            # A tag for the output files. 

# Concatenates multiple runs into one big data file for PCA

pdir="projects/neuromat/00-DATA/2013-05-23-ghislain"

if [[ ${filter} -gt 0 ]]; then
  subdir="flt-runs" # Source subdir of ${pdir}.
  stag="_f"         # Tag of source files.
  ftag="_y"         # Tag of destination files.
else
  subdir="raw-runs" # Source subdir of ${pdir}.
  stag=""           # Tag of source files.
  ftag="_x"         # Tag of destination files.
fi

afile="`echo ${pdir}/${subdir}/s${subject}_r001${stag}.txt`"
echo "reading header data from ${afile}"  1>&2

# Get channel names and count:
channels=( `cat ${afile} | gawk '/^channels = /{ $2=""; $1=""; print; exit(0); }'`  )
nc=${#channels[@]}
if [[ ${nc} -lt 1 ]]; then
   echo "** bad {nc} in \"${afile}\"" 1>&2 ; exit 1
fi

# Get electrode count {ne}:
ne=( `cat ${afile} | gawk '/^ne = /{ print $3; exit(0); }'`  )
if [[ ${#ne[@]} -ne 1 ]]; then
   echo "** bad {ne} in \"${afile}\"" 1>&2 ; exit 1
fi
ne="${ne[0]}"

# Get sampling frequency:
fsmp=( `cat ${afile} | gawk '/^fsmp = /{ print $3; exit(0); }'`  )
if [[ ${#fsmp[@]} -ne 1 ]]; then
   echo "** bad {fsmp} in \"${afile}\"" 1>&2 ; exit 1
fi
fsmp="${fsmp[0]}"

ofile="data/s${subject}_${otag}${ftag}.txt"
echo "writing multi-run dataset to ${ofile}" 1>&2

# Write a minimal header of output file
rm -fv ${ofile}
echo "nc = ${nc}" >> ${ofile}
echo "ne = ${ne}" >> ${ofile}
echo "channels = ${channels[*]}" >> ${ofile}
echo "fsmp = ${fsmp}" >> ${ofile}

# Concatenate the runs with apodizing mask:
touch .DUMMY
run=${rini};
while [[ ${run} -le ${rfin} ]]; do
  xrun=`printf "%03d" "${run}"` 
  rfile="${pdir}/${subdir}/s${subject}_r${xrun}${stag}.txt"
  
  echo "  ${rfile}" 1>&2
  # Append ${rfile}: 
  cat ${rfile} \
    | grep -v -e '=' \
    | soften_dataset_ends.gawk \
        -v nc="${nc}" -v ne="${ne}" -v fsmp="${fsmp}" \
    >> ${ofile}
  
  run=$(( ${run} + 1 ))
done
