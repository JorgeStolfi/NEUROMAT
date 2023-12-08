# Last edited on 2023-11-02 15:26:08 by stolfi

PROG := nmeeg_plot_electrode_histogram.sh
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# Format: {DATASET}:{CHNAME}

TESTS := \
  synth_s:C002 \
  s003_r009_e20_r:F7

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
CHNAME :=  ${word 2,${TEST_ARGS}}

SHOW := SHOW

HIST_DATA_FILE := data/${DATASET}_hist_${CHNAME}
HIST_PLOT_FILE := out/${DATASET}_hist_${CHNAME}

do-single-test: ${PROGDIR}/${PROG}
	rm -f out/${DATASET}*
	${PROGDIR}/${PROG} \
          ${SHOW} \
          ${CHNAME} \
          ${HIST_DATA_FILE} \
          ${HIST_PLOT_FILE}

# End ${TEST} section.
# ======================================================================
endif
