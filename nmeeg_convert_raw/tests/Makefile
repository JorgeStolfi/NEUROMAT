# Last edited on 2023-11-29 11:24:59 by stolfi

PROG := nmeeg_convert_raw
PROGDIR = ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

SHOW := SHOW

TESTS := 13 14
 
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
SUBJECT :=   ${word 1,${TEST_ARGS}}

DATASET_NAME := S${SUBJECT}

DATA_DIR := projects/neuromat/00-DATA/2013-09-04-ghislain-128/raw

RAW_EEG_FILE := ${DATA_DIR}/${DATASET_NAME}.raw

OUT_DIR := out
OUT_PREFIX := ${OUT_DIR}/${DATASET_NAME}-chunk

OUT_FILE := ${OUT_PREFIX}.txt

# Max voltage (uV) for plot scales
VMAX_FULL := 15000
VMAX_EVTS := 5

do-single-test: ${PROGDIR}/${PROG} ${RAW_EEG_FILE}
	rm -f ${OUT_PREFIX}.txt
	cat ${RAW_EEG_FILE} \
          | ${PROGDIR}/${PROG} \
              -unit 1.0 \
              -skipFrames 550000 \
              -copyFrames   5000 \
              -sourceFile ${DATASET_NAME}.raw \
              -subject ${SUBJECT} \
          > ${OUT_FILE}
	@echo "== plotting the input data =="
	nmeeg_plot_channels.sh   ${SHOW} ${OUT_PREFIX} 0 ${VMAX_FULL}   0 9999     0 0  1400 500
	nmeeg_plot_channels.sh   ${SHOW} ${OUT_PREFIX} 0 ${VMAX_EVTS} 129 9999     0 0  1400 500

# ======================================================================
endif
