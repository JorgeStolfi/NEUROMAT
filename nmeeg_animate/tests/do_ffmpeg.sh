#! /bin/bash
# Last edited on 2021-08-29 00:48:55 by stolfi

frameDir="$1"; shift
outFile="$1"; shift

ffmpeg \
  -framerate 30 \
  -start_number 0 \
  -i "${frameDir}/f%06d.png" \
  -y \
  -r 30 \
  -vcodec libx264 \
  -crf 25 \
  -pix_fmt yuv420p \
  ${outFile}
  
#  -vf "pad=width=540:height=720:x=45:y=22:color=black" \
