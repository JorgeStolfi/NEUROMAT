# Last edited on 2023-12-05 10:25:00 by stolfi

PROG := nmeeg_plot_channels.sh
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# Format: {DATASET}:{N_SKIP}:{N_USE}:{FSMP}
 
TESTS := \
  s001_r005:000000:004800:600 \
  s001_r005:001510:000580:600

ALL_TESTS := \
  s001_r005:000000:004800:600 \
  s001_r005:001510:000580:600

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
N_SKIP :=   ${word 2,${TEST_ARGS}}
N_USE :=    ${word 3,${TEST_ARGS}}
FSMP :=     ${word 4,${TEST_ARGS}}

VMAX := 100

SHOW := SHOW

EEG_FILE := data/${DATASET}.txt

OUT_PREFIX := out/${DATASET}_${N_SKIP}_${N_USE}

do-single-test: ${PROGDIR}/${PROG}
	rm -f ${OUT_PREFIX}*
	${PROGDIR}/${PROG} \
          ${SHOW} \
          data/${DATASET} 0 ${VMAX} \
          5 7 \
          ${N_SKIP} ${N_USE} \
          1400 500 \
          -marker 20 0.02 -80

# ======================================================================
endif
