# Last edited on 2014-01-23 03:03:17 by stolfilocal

PROG := nmeeg_animate
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# Format: {DATASET}:{NE}:{N_SKIP}:{N_USE}:{RSTEP}:{B_TYPE}
 
TESTS := \
  s013_r00112_f:128:000000:000700:1:03:256:256

SINGLE_TESTS := \
  s013_r00112_f:128:000100:000001:1:00:256:256 \
  s013_r00112_f:128:000100:000001:1:01:256:256 \
  s013_r00112_f:128:000100:000001:1:02:256:256

QUICK_TESTS := \
  s013_r00112_f:128:000100:000150:1:00:256:256 \
  s013_r00112_f:128:000100:000150:1:01:256:256 \
  s013_r00112_f:128:000100:000150:1:02:256:256

ALL_TESTS := \
  s001_r005_r:20:000500:003800:5:02:280:320 \
  s001_r005_r:20:000500:003800:5:01:280:320 \
  s001_r005_r:20:001000:001800:1:01:280:320 \
  \
  s013_r00112_f:128:000000:000000:1:01:256:256 \
  s013_r00112_f:128:000000:000000:1:02:256:256

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
NE :=       ${word 2,${TEST_ARGS}}
N_SKIP :=   ${word 3,${TEST_ARGS}}
N_USE :=    ${word 4,${TEST_ARGS}}
RSTEP :=    ${word 5,${TEST_ARGS}}
B_TYPE :=   ${word 6,${TEST_ARGS}}
NX :=       ${word 7,${TEST_ARGS}}
NY :=       ${word 8,${TEST_ARGS}}

SHOW := SHOW

EEG_FILE := data/${DATASET}.txt

B_PREFIX := data/${NX}_${NY}_e${NE}

OUT_PREFIX := out/${DATASET}_${N_SKIP}_${N_USE}_s${RSTEP}_b${B_TYPE}

GIF_FILE := ${OUT_PREFIX}.gif

do-single-test: ${PROGDIR}/${PROG} data/${DATASET}.txt
	rm -f out/*_f??????.png ${OUT_PREFIX}_* ${GIF_FILE}
	${PROGDIR}/${PROG} \
          -basis ${B_PREFIX} ${B_TYPE} \
          -skip ${N_SKIP} -read ${N_USE} -step ${RSTEP} \
          -outPrefix ${OUT_PREFIX} \
          < ${EEG_FILE}
	convert \
          ${B_PREFIX}_msk.png ${B_PREFIX}_msk.png ${B_PREFIX}_msk.png \
          `ls ${OUT_PREFIX}_f??????.png | sort` -loop 2 -delay 5/100 \
          ${GIF_FILE}
	rm -f ${OUT_PREFIX}_f??????.png
	if [[ "${SHOW}" == "SHOW" ]]; then nomacs ${GIF_FILE}; fi
        
# ----------------------------------------------------------------------
# Obtaining the test data:

# Data source for 20-electrode raw runs:
RAW_E20_RUNS_SOURCE_DIR := ../../nmeeg_split_e20_gh2012/tests/out/

data/%_r.txt: ${RAW_E20_RUNS_SOURCE_DIR}/%.txt
	cp -av ${RAW_E20_RUNS_SOURCE_DIR}/$*.txt data/$*_r.txt

# Data source for 20-electrode filtered runs:
FLT_E20_RUNS_SOURCE_DIR := ../../nmeeg_filter/tests/out/

data/%_f.txt: ${FLT_E20_RUNS_SOURCE_DIR}/%_f.txt
	cp -av ${FLT_E20_RUNS_SOURCE_DIR}/$*_f.txt data/$*_f.txt

# Data source for 128-electrode filtered runs:
FLT_E128_RUNS_SOURCE_DIR := projects/neuromat/00-DATA/2013-09-04-ghislain-128/flt-runs-B

data/s013_%_f.txt: ${FLT_E128_RUNS_SOURCE_DIR}/s013_%.txt
	cp -av ${FLT_E128_RUNS_SOURCE_DIR}/s013_$*.txt data/s013_$*_f.txt


# ======================================================================
endif
