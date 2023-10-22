# Last edited on 2021-08-30 18:20:51 by stolfi

PROG := nmeeg_make_noise_data.gawk
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# Format: {NT}
 
TESTS := \
  3500 \

ALL_TESTS :=  \
  3500

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
NT :=    ${word 1,${TEST_ARGS}}

OUT_PREFIX := out/test_${TYPE}_${NT}_${FSMP}
OUT_TXT_FILE := ${OUT_PREFIX}.txt
OUT_PNG_FILE := ${OUT_PREFIX}.png

do-single-test: ${PROGDIR}/${PROG}
	rm -fv ${OUT_TXT_FILE}
	${PROGDIR}/${PROG} \
            -v nt=${NT} \
          > ${OUT_TXT_FILE}
	nmeeg_plot_channels.sh \
          SHOW \
          ${OUT_PREFIX} 0 \
          41 0 -1 \
          0 0 \
          500

# ======================================================================
endif
