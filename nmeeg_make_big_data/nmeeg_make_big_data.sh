#! /bin/bash
# Last edited on 2023-11-04 09:20:07 by stolfi

# Creates an EEG dataset test file by extracting a largish chunk of
# one of Ghilslain's datasets, or filtering such a chunk.

# Arguments:
subject=$(echo "$1" | sed -e 's:^0*\([0-9]\):\1:'); shift  # Single-digit subject ID.
block="$1"; shift                                          # Block id ("1", "2", or "1-2").
filter=$(echo "$1" | sed -e 's:^0*\([0-9]\):\1:'); shift   # Use filtered (1) or unfiltered (0) runs?
otag="$1"; shift                                           # Output tag.
nskip=$(echo "$1" | sed -e 's:^0*\([0-9]\):\1:'); shift;   # Index of initial frame (from 0).
nread=$(echo "$1" | sed -e 's:^0*\([0-9]\):\1:'); shift;   # Number of frames to extract.

# If ${filter} is 0,extracts from a raw dataset 
# "s${subject}_bl${block}.txt" a chunk 
# consecutive frames, disregarding runs etc.
# Writes "data/s00${subject}_${otag}_x.txt"
# If ${filter} is 1, filters the previously extracted dataset
# yielding "data/s00${subject}_${otag}_y.txt"

# Variables:
pdir="projects/neuromat/00-DATA/2013-05-23-ghislain"

rawfile="data/s00${subject}_${otag}_x.txt"

channels=( `cat ${pdir}/channel_names.txt | gawk '/^[ ]*[0-9]/{ print $3; }'` )

if [[ ${filter} -eq 0 ]]; then 

  # Extract frame subset and write to ${rawfile}:
  
  rm -fv ${rawfile}

  echo "nc = 21" >> ${rawfile}
  echo "ne = 20" >> ${rawfile}
  echo "channels = ${channels[*]}" >> ${rawfile}
  echo "fsmp = 600" >> ${rawfile}

  cat ${pdir}/raw/s${subject}_bl${block}*.txt.gz \
    | gunzip \
    | tr -d '\015' \
    | head -n $(( ${nskip} + ${nread} )) \
    | tail -n $(( ${nread} )) \
    >> ${rawfile}

else
  # Filter the previously extracted ${rawfile}
  filfile="data/s00${subject}_${otag}_y.txt"
  filpref=${filfile%%.*}
  
  # Filter parameters (withot invert flag):
  filparms=( "G" 0.05 0.15 25.0 30.0 )

  cat ${rawfile} \
    | nmeeg_filter \
        ${filpref} 0 0  \
        ${filparms[@]} 0

  mv -v ${filpref}_f.txt ${filfile}
fi
