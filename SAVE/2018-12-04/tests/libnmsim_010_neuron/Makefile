# Last edited on 2017-07-21 14:55:42 by stolfilocal

PROG := nmsim_test_010_neuron
 
TEST_LIB := libnmsim.a

TEST_LIB_DIR := ../..

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
all: do-test

do-test: ${PROG}
	${PROG}
