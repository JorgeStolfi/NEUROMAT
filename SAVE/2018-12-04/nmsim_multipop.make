# Last edited on 2016-12-14 00:23:01 by stolfilocal

PROG := nmsim_multipop

JS_LIBS := \
  libnmsim.a \
  libjs.a

USE_GL := NO
USE_X11 := NO
  
IGNORE :=

all: build install

include GENERIC-PROGS.make

