# Last edited on 2013-06-03 23:11:37 by stolfilocal

PROG = nmeeg_animate

JS_LIBS := \
  libneuro.a \
  libimg.a \
  libgeo.a \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib64/libjpeg.so \
  /usr/lib64/libpng.a  \
  /usr/lib64/libz.so

include ${STOLFIHOME}/programs/c/GENERIC-PROGS.make
