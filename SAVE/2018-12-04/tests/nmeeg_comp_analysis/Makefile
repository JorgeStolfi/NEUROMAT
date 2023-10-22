# Last edited on 2013-12-04 11:13:16 by stolfilocal

PROG := nmeeg_comp_analysis
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# Format: {SUBJECT}:{RUN}:{SCRTAG}
 
TESTS := \
  013:00229:blink \
  014:00425:blink \
  013:00229:blinkosc \
  014:00425:blinkosc
 
BLINK_TESTS := \
  013:00229:blink \
  014:00425:blink \
  013:00229:blinkosc \
  014:00425:blinkosc

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
RUNID :=      ${word 2,${TEST_ARGS}}
SCRTAG :=     ${word 3,${TEST_ARGS}}

DATA_DIR := data
OUT_DIR := out

DATASET := s${SUBJID}_r${RUNID}
EEG_FILE := ${DATA_DIR}/${DATASET}.txt
SCRIPT := do_test_${SCRTAG}.sh
OUT_PREFIX := ${OUT_DIR}/${DATASET}_${SCRTAG}

do-single-test: ${PROGDIR}/${PROG} ${EEG_FILE}
	rm -f ${OUT_PREFIX}{,_*}.{png,txt}
	./${SCRIPT} ${SUBJID} ${RUNID} ${EEG_FILE} ${OUT_PREFIX}
	if [[ "/${SHOW}" == "/SHOW" ]]; then \
          echo "plotting the decomposed signals over time" 1>&2 ; \
          for ff in ${OUT_PREFIX}_*.txt ; do \
            prefix=$${ff%.txt} ; \
            nmeeg_plot_channels.sh SHOW $${prefix} 0 180 0 -1 0 0 ; \
          done ; \
          for ff in ${OUT_PREFIX}.txt ; do \
            prefix=$${ff%.txt} ; \
            nmeeg_plot_channels.sh SHOW $${prefix} 0 50 0 -1 0 0 ; \
          done ; \
        fi

later:
	echo "computing spectra of principal components" 1>&2 ; \
          nmeeg_spectrum ${OUT_PREFIX}_pca 0 0 < ${OUT_PREFIX}_pca.txt ; \
          nmeeg_plot_spectra.sh ${SHOW} ${OUT_PREFIX}_pca_pwr 100 0 19 480 'linespoints pt 7'
	if [[ "/${ANIMATE}" == "/ANIMATE" ]]; then \
          echo "animating the principal component patterns" 1>&2 ; \
          animate_components.sh ${SHOW} ${OUT_PREFIX} ${N_USE} ${RSTEP} ; \
        fi
        
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# Obtaining the test data:

# ----------------------------------------------------------------------
# Data source for raw runs:
RAW_RUNS_SOURCE_DIR := ../../nmeeg_split_e20_gh2012/tests/out/

${DATA_DIR}/%_r.txt: ${RAW_RUNS_SOURCE_DIR}/%.txt
	cp -av ${RAW_RUNS_SOURCE_DIR}/$*.txt ${DATA_DIR}/$*_r.txt

# ----------------------------------------------------------------------
# Data source for filtered runs:
FILT_RUNS_SOURCE_DIR := ../../nmeeg_filter/tests/out/

${DATA_DIR}/%_f.txt: ${FILT_RUNS_SOURCE_DIR}/%_r_f.txt
	cp -av ${FILT_RUNS_SOURCE_DIR}/$*_r_f.txt ${DATA_DIR}/$*_f.txt

# ----------------------------------------------------------------------
# Raw data spanning multiple runs:
PROJECT_DIR := projects/neuromat/00-DATA/2013-05-23-ghislain

${DATA_DIR}/s004_longA_x.txt: ${PROJECT_DIR}/raw/s4_bl1.txt.gz ~/bin/nmeeg_make_big_data.sh
	nmeeg_make_big_data.sh 4 1 0 longA 30000 18000

${DATA_DIR}/s004_longA_y.txt: ${DATA_DIR}/s004_longA_x.txt ~/bin/nmeeg_make_big_data.sh
	nmeeg_make_big_data.sh 4 1 1 longA 30000 18000

${DATA_DIR}/s006_longA_x.txt: ${PROJECT_DIR}/raw/s6_bl1.txt.gz ~/bin/nmeeg_make_big_data.sh
	nmeeg_make_big_data.sh 6 1 0 longA 10000 18000

${DATA_DIR}/s006_longA_y.txt: ${DATA_DIR}/s006_longA_x.txt ~/bin/nmeeg_make_big_data.sh
	nmeeg_make_big_data.sh 6 1 1 longA 10000 18000

# ----------------------------------------------------------------------
# Concatenation of many runs:
${DATA_DIR}/s004_mrunA_x.txt: ${PROJECT_DIR}/raw-runs/s004_r001.txt ./make_mrun_data.sh
	./make_mrun_data.sh 004 0 10 40 mrunA

${DATA_DIR}/s004_mrunA_y.txt: ${PROJECT_DIR}/flt-runs/s004_r001_f.txt ./make_mrun_data.sh
	./make_mrun_data.sh 004 1 10 40 mrunA

# ======================================================================
endif
	
