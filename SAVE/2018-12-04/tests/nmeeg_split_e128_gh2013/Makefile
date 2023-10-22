# Last edited on 2013-11-13 22:02:45 by stolfilocal

PROG := nmeeg_split_e128_gh2013
PROGDIR = ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# ${TESTS} format: {DATASET}:{SUBJECT}:{FIRST_RUN}:{IN_FRAMES}:{SAMPLE_TAG}
# where {SUBJECT} is 3 digits decimal.

TESTS := \
  S14:014:01:30000:r00103

ALL_TESTS := \
  S13:013:01:30000:r00101 \
  S14:014:01:30000:r00101
 
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
DATASET_NAME :=   ${word 1,${TEST_ARGS}}
SUBJECT :=        ${word 2,${TEST_ARGS}}
FIRST_RUN :=      ${word 3,${TEST_ARGS}}
IN_FRAMES :=      ${word 4,${TEST_ARGS}}
SAMPLE_TAG :=     ${word 5,${TEST_ARGS}}

RAW_EEG_FILE := data/${DATASET_NAME}.raw
FSMP := 500

OUT_DIR := out
OUT_PREFIX := ${OUT_DIR}/s${SUBJECT}

TMP_FILE := /tmp/.eeg.txt

FRAMES_PRE := 1800
FRAMES_POS := 1800
FRAMES_PER_RUN := ${shell echo '${FRAMES_PRE} + ${FRAMES_POS}' | bc -lq}

SAMPLE_NAME = ${OUT_PREFIX}_${SAMPLE_TAG}
SAMPLE_FILE = ${SAMPLE_NAME}.txt}

do-single-test: ${PROGDIR}/${PROG} ${RAW_EEG_FILE} ${RTY_FILE}
	rm -f ${OUT_PREFIX}_r?????.txt
	cat ${RAW_EEG_FILE} \
          | nmeeg_convert_raw \
             -unit 1.00 \
             -skipFrames 0 \
             -copyFrames ${IN_FRAMES} \
             -sourceFile ${DATASET_NAME} \
             -subject ${SUBJECT} \
            > ${TMP_FILE}
	${PROGDIR}/${PROG} \
              -firstBlock 01 \
              -firstRun 01 \
              -runsPerBlock 32 \
              -framesPerRun ${FRAMES_PRE} ${FRAMES_POS} \
              -stimulusChannel HB   4.00  50.000 HB \
              -stimulusChannel HS   4.00  50.000 HS \
              -stimulusChannel QB   4.00  50.000 QB \
              -stimulusChannel QS   4.00  50.000 QS \
              -fixationChannel fix1 4.00  50.000 \
              -outDir ${OUT_DIR} \
            < ${TMP_FILE}
	wc -l ${OUT_PREFIX}_r?????.txt
	nmeeg_plot_channels.sh SHOW ${SAMPLE_NAME} ${FSMP} 100 0 -1 0 ${FRAMES_PER_RUN} 

# ----------------------------------------------------------------------
# Obtaining the test data:

# Data source beginning with "projects"
DATA_SOURCE_DIR := projects/neuromat/00-DATA/2013-09-04-ghislain-128

data/%.raw: ${DATA_SOURCE_DIR}/raw/%.raw
	( cd data && ln -s ${DATA_SOURCE_DIR}/raw/$*.raw )

# ======================================================================
endif
