# Last edited on 2023-12-05 10:29:29 by stolfi

PROG := nmeeg_plot_bis_channels.sh
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# Format: {DATASET}:{N_SKIP}:{N_USE}:{FSMP}
 
TESTS := \
  s001_r005:000000:004800:600 \
  s001_r005:001510:001780:600

ALL_TESTS := \
  s001_r005:000000:004800:600 \
  s001_r005:001510:001780:600

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

N_SKIP0 := ${shell echo "${N_SKIP} + 250" | bc -lq}
N_USE0 := 90

N_SKIP1 := ${N_SKIP}
N_USE1 := ${N_USE}

VMAX := 100

SHOW := SHOW

EEG_FILE := data/${DATASET}.txt

OUT_PREFIX := out/${DATASET}_${N_SKIP}_${N_USE}

do-single-test: ${PROGDIR}/${PROG}
	rm -f ${OUT_PREFIX}*.png
	${PROGDIR}/${PROG} \
          ${SHOW} \
          data/${DATASET} ${FSMP} ${VMAX} \
          5 7 \
          ${N_SKIP0} ${N_USE0} 400 \
          ${N_SKIP1} ${N_USE1} 1200 \
          500 \
          -outPrefix ${OUT_PREFIX} \
          -marker 20 0.02 -80

# ======================================================================
endif
