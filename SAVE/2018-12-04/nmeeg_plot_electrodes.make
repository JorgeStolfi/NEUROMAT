# Last edited on 2013-09-30 23:47:28 by stolfilocal

PROG = nmeeg_plot_electrodes

JS_LIBS := \
  libneuro.a \
  libgeo.a \
  libps.a \
  libjs.a \
  

include ${STOLFIHOME}/programs/c/GENERIC-PROGS.make
