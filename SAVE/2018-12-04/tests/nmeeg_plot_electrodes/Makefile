# Last edited on 2013-10-01 00:23:24 by stolfilocal

PROG := nmeeg_plot_electrodes
PROGDIR = ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

TESTS := 020 128

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
NELEC :=   ${word 1,${TEST_ARGS}}

RADIUS := 100

OUT_DIR := out
OUT_PREFIX := ${OUT_DIR}/fig-${NELEC}

do-single-test: ${PROGDIR}/${PROG}
	rm -f ${OUT_PREFIX}-els.eps
	${PROGDIR}/${PROG} \
            -radius ${RADIUS} \
            -nElec ${NELEC} \
            -outName ${OUT_PREFIX}
	gv ${OUT_PREFIX}-els.eps

# ======================================================================
endif
