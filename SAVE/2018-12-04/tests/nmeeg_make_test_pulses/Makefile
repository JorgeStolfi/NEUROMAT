# Last edited on 2013-11-16 03:17:24 by stolfilocal

PROG := nmeeg_make_test_pulses.gawk
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# Format: {DATASET}:{N_USE}:{FSMP}:{VARIED}
 
TESTS := \
  pulsesB:3600:600:1 \
  pulsesA:3600:600:0

ALL_TESTS := \
  pulsesA:3600:600:0 \
  pulsesB:3600:600:1

.PHONY:: do-all-tests do-single-test

all: do-all-tests

do-all-tests:
	for tt in ${TESTS} ; do \
          ${MAKE} TEST=$$tt do-single-test ; \
        done

ifneq "/${TEST}" "/"
# ======================================================================
# Recursive section, client must define ${TEST}
  
TEST_ARGS := ${subst :, ,${TEST}}
DATASET :=  ${word 1,${TEST_ARGS}}
N_USE :=    ${word 2,${TEST_ARGS}}
FSMP :=     ${word 3,${TEST_ARGS}}
VARIED :=   ${word 4,${TEST_ARGS}}

OUT_PREFIX := out/${DATASET}

do-single-test: ${PROGDIR}/${PROG}
	rm -f ${OUT_PREFIX}_*
	${PROGDIR}/${PROG} \
            -v varied=${VARIED} \
            -v nt=${N_USE} \
          > ${OUT_PREFIX}.txt
	nmeeg_plot_channels.sh \
          SHOW \
          ${OUT_PREFIX} 0 \
          160 0 -1 \
          0 0

# ======================================================================
endif
