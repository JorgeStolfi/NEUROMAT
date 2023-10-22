# Last edited on 2013-12-06 03:42:36 by stolfilocal

PROG = nmeeg_make_imgbasis

JS_LIBS := \
  libneuro.a \
  libimg.a \
  libgeo.a \
  libjs.a
  
OTHER_LIBS :=

include ${STOLFIHOME}/programs/c/GENERIC-PROGS.make
