# Last edited on 2013-06-03 03:22:20 by stolfilocal

PROG = nmeeg_filter

JS_LIBS := \
  libneuro.a \
  libgeo.a \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib64/libfftw3.so \

include ${STOLFIHOME}/programs/c/GENERIC-PROGS.make
