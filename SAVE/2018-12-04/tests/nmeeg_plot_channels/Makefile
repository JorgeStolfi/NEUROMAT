# Last edited on 2013-06-12 11:25:01 by stolfilocal

PROG := nmeeg_plot_channels.sh
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# Format: {DATASET}:{N_SKIP}:{N_USE}:{FSMP}
 
TESTS := \
  s001_r005:000000:004800:600 \
  s001_r005:003000:000600:600

ALL_TESTS := \
  s001_r005:000000:004800:600 \
  s001_r005:003000:000600:600

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

SHOW := SHOW

EEG_FILE := data/${DATASET}.txt

OUT_PREFIX := out/${DATASET}_${N_SKIP}_${N_USE}

do-single-test: ${PROGDIR}/${PROG}
	rm -f ${OUT_PREFIX}_*
	${PROGDIR}/${PROG} \
          ${SHOW} \
          data/${DATASET} ${FSMP} \
          160 0 20 \
          ${N_SKIP} ${N_USE}

# ======================================================================
endif
