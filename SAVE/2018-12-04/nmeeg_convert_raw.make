# Last edited on 2013-10-01 15:57:37 by stolfilocal

PROG = nmeeg_convert_raw

JS_LIBS := \
  libneuro.a \
  libjs.a
  
OTHER_LIBS :=

include ${STOLFIHOME}/programs/c/GENERIC-PROGS.make
