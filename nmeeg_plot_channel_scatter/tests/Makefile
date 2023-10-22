# Last edited on 2021-09-02 20:36:30 by stolfi

PROG := nmeeg_plot_channel_scatter.sh
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# Format: {DATASET}:{STEP}:{XNAME}:{YNAME}
 
TESTS := \
  r00624_AMP:002:P000:P003 \
  r00624_AMP:002:P004:P003

ALL_TESTS := \
  r00624_AMP:002:P000:P003 \
  r00624_AMP:002:P000:P003

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
STEP :=     ${word 2,${TEST_ARGS}}
XNAME :=    ${word 3,${TEST_ARGS}}
YNAME :=    ${word 4,${TEST_ARGS}}

SHOW := SHOW

EEG_FILE := data/${DATASET}.txt

OUT_FILE := out/${DATASET}_${XNAME}_${YNAME}_scatter.png

do-single-test: ${PROGDIR}/${PROG}
	rm -f ${OUT_FILE}
	${PROGDIR}/${PROG} \
          ${SHOW} \
          data/${DATASET}.txt \
          ${STEP} \
          ${XNAME} 150 \
          ${YNAME} 150

# ======================================================================
endif
