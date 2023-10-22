# Last edited on 2017-10-17 21:25:10 by jstolfi

PROGDIR := ..
PROG := nmeeg_cleanup

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# Format: {SUBJECT}:{RUNS}

TESTS := \
  013:00101,00102,00305

.PHONY:: do-all-tests do-single-test get-data

SHOW := SHOW
ANIMATE := NOANIMATE

all: do-all-tests

do-all-tests:
	for tt in ${TESTS} ; do \
          ${MAKE} TEST=$$tt do-single-test do-single-show ; \
        done

ifneq "/${TEST}" "/"
# ======================================================================
# Recursive section, client must define ${TEST}
  
TEST_ARGS := ${subst :, ,${TEST}}
SUBJID :=     ${word 1,${TEST_ARGS}}
RUNS :=       ${word 2,${TEST_ARGS}}

DATA_DIR := data
OUT_DIR := out

OUT_PREFIX := ${OUT_DIR}/s${SUBJID}
IN_PREFIX := ${DATA_DIR}/s${SUBJID}
COMMA := ,

DATA_FILES := ${addsuffix .txt, ${addprefix ${IN_PREFIX}_r, ${subst ${COMMA}, ,${RUNS}}}}

SCRIPT := do_test.sh

do-single-test: ${PROGDIR}/${PROG} ${DATA_FILES}
	rm -f ${OUT_PREFIX}_r*.{png,txt}
	./${SCRIPT} ${IN_PREFIX} "${RUNS}" ${OUT_PREFIX} ${SHOW} ${ANIMATE}
        
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# Obtaining the test data:

# ----------------------------------------------------------------------
# Data source for filtered runs:
FLT_RUNS_SOURCE_DIR := ~/projects/neuromat/00-DATA/2013-09-04-ghislain-128/flt-runs-B

${DATA_DIR}/%.txt: ${FLT_RUNS_SOURCE_DIR}/%.txt
	cp -av ${FLT_RUNS_SOURCE_DIR}/$*.txt ${DATA_DIR}/

# ======================================================================
endif
	

