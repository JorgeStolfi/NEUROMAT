# Last edited on 2023-11-29 11:21:11 by stolfi

PROG := nmeeg_split_e20_gh2012
PROGDIR = ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# ${TESTS} format: {DATASET}:{SUBJECT}:{FIRST_RUN}:{FSMP}:{SAMPLE_TAG}
# where {SUBJECT} is 3 digits decimal.

TESTS := \
  s1_bl1:001:001:600:r005 \
  s2_bl1:002:001:600:r005 \
  s5_bl1-2:005:001:600:r090 \
  s8_bl1:008:001:600:r005

ALL_TESTS := \
  s1_bl1:001:001:600:r005 \
  s2_bl1:002:001:600:r005 \
  s5_bl1-2:005:001:600:r090 \
  s8_bl1:008:001:600:r005
 
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
FSMP :=           ${word 4,${TEST_ARGS}}
SAMPLE_TAG :=     ${word 5,${TEST_ARGS}}

RAW_EEG_FILE := data/${DATASET_NAME}.txt.gz

RTY_FILE := data/${DATASET_NAME}.rty

CHANNELS := F7 T3 T5 Fp1 F3 C3 P3 O1 F8 T4 T6 Fp2 F4 C4 P4 O2 Fz Cz Pz Oz TR

OUT_DIR := out
OUT_PREFIX := ${OUT_DIR}/s${SUBJECT}

PULSES_PER_RUN := 4
FRAMES_PER_RUN := 6600

SAMPLE_NAME = ${OUT_PREFIX}_${SAMPLE_TAG}
SAMPLE_FILE = ${SAMPLE_NAME}.txt}

do-single-test: ${PROGDIR}/${PROG} ${RAW_EEG_FILE} ${RTY_FILE}
	rm -f ${OUT_PREFIX}_r???.txt
	runtypes=( `cat ${RTY_FILE} | sed -e 's:nonBio:N:g' -e 's:Bio:B:g'` ) ; \
	cat ${RAW_EEG_FILE} | gunzip \
          | ${PROGDIR}/${PROG} \
              -source ${DATASET_NAME} \
              -subject ${SUBJECT} \
              -outDir ${OUT_DIR} \
              -framesPerRun ${FRAMES_PER_RUN} \
              -firstRun ${FIRST_RUN} \
              -pulsesPerRun ${PULSES_PER_RUN} \
              -channels ${CHANNELS} \
              -trigger TR -nElectrodes 20 \
              -fSampling ${FSMP} \
              -runTypes $${runtypes[@]} \
              -marker "SOn" 0 1  -12  600 30.0 \
              -marker "Occ" 1 1 -612 1800 50.0 \
              -marker "RAp" 2 1 -312  900 60.0
	wc -l ${OUT_PREFIX}_r???.txt
	nmeeg_plot_channels.sh SHOW ${SAMPLE_NAME} ${FSMP} 110  0 9999  0 ${FRAMES_PER_RUN}  1400 500

# ----------------------------------------------------------------------
# Obtaining the test data:

# Data source beginning with "projects"
DATA_SOURCE_DIR := projects/neuromat/00-DATA/2013-05-23-ghislain

data/%.txt.gz: ${DATA_SOURCE_DIR}/raw/%.txt.gz
	( cd data && ln -s ${DATA_SOURCE_DIR}/raw/$*.txt.gz )

data/%.rty: ${DATA_SOURCE_DIR}/raw/%_seq.log ${DATA_SOURCE_DIR}/runtypes_from_log.gawk
	( cd data && cat ${DATA_SOURCE_DIR}/raw/$*_seq.log \
          | ${DATA_SOURCE_DIR}/runtypes_from_log.gawk \
          > $*.rty )

# ======================================================================
endif
