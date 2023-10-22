# Last edited on 2013-11-13 21:05:32 by stolfilocal

PROG = nmeeg_spectrum

JS_LIBS := \
  libneuro.a \
  libgeo.a \
  libjs.a
  
OTHER_LIBS := \
  /usr/lib64/libfftw3.a

include ${STOLFIHOME}/programs/c/GENERIC-PROGS.make
