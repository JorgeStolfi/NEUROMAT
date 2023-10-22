# Last edited on 2013-12-02 04:24:41 by stolfilocal

PROG = nmeeg_correl

JS_LIBS := \
  libneuro.a \
  libgeo.a \
  libjs.a
  
OTHER_LIBS :=

include ${STOLFIHOME}/programs/c/GENERIC-PROGS.make
