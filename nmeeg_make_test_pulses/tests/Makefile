# Last edited on 2021-08-26 18:59:42 by stolfi

PROG := nmeeg_make_test_pulses.gawk
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# Format: {TYPE}:{NT}:{FSMP}
 
TESTS := \
  3:3500:500 \
  2:3500:500 \
  1:3500:500 \
  0:3500:500 \

ALL_TESTS := \
  0:3500:500 \
  1:3500:500 \
  2:3500:500 \
  3:3500:500 

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
TYPE :=  ${word 1,${TEST_ARGS}}
NT :=    ${word 2,${TEST_ARGS}}
FSMP :=  ${word 3,${TEST_ARGS}}

OUT_PREFIX := out/test_${TYPE}_${NT}_${FSMP}
OUT_TXT_FILE := ${OUT_PREFIX}.txt
OUT_PNG_FILE := ${OUT_PREFIX}.png

do-single-test: ${PROGDIR}/${PROG}
	rm -fv ${OUT_TXT_FILE}
	${PROGDIR}/${PROG} \
            -v type=${TYPE} \
            -v nt=${NT} \
            -v fsmp=${FSMP} \
          > ${OUT_TXT_FILE}
	nmeeg_plot_channels.sh \
          SHOW \
          ${OUT_PREFIX} 0 \
          160 0 -1 \
          0 0 \
          500

# ======================================================================
endif
