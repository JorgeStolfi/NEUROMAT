# Last edited on 2023-11-29 11:22:01 by stolfi

PROG := nmeeg_spectrum
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# Format: {DATASET}
 
TESTS := \
  s013_r00101_e128_r

ALL_TESTS := \
  s013_r00101_e128_r \
   \
  s005_r044_e20_r \
  s005_r044_e20_f \
   \
  s008_r015_e20_r \
  s008_r015_e20_f \
   \
  s008_r049_e20_r \
  s008_r049_e20_f \
   \
  syntA_s \
   \
  pulseA_s

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

SHOW := SHOW

EEG_FILE := data/${DATASET}.txt

NT :=    ${shell egrep -e '^ *nt = ' ${EEG_FILE} | gawk '//{ print $$3; }' | head -n 1}
NC :=    ${shell egrep -e '^ *nc = ' ${EEG_FILE} | gawk '//{ print $$3; }' | head -n 1}
NE :=    ${shell egrep -e '^ *ne = ' ${EEG_FILE} | gawk '//{ print $$3; }' | head -n 1}
FSMP :=  ${shell egrep -e '^ *fsmp = ' ${EEG_FILE} | gawk '//{ print $$3; }' | head -n 1 | sed -e 's:[.][0]*$$::g'}

FMAX := ${shell echo 'scale=0; ${FSMP}/2' | bc -lq}

OUT_PREFIX := out/${DATASET}

do-single-test: ${PROGDIR}/${PROG} ${EEG_FILE}
	rm -f ${OUT_PREFIX}_* ${OUT_PREFIX}.*
	nmeeg_plot_channels.sh   ${SHOW} data/${DATASET} ${FSMP} 420  0 $$(( ${NC} - 1 ))  0 ${NT}  1400 500
	${PROGDIR}/${PROG} \
          ${OUT_PREFIX} 0 ${NT} \
          < ${EEG_FILE}
	nmeeg_plot_spectra.sh    ${SHOW} ${OUT_PREFIX}_pwr ${FMAX} 0 $$(( ${NE} - 1 )) 500 'linespoints pt 7'
        
# ----------------------------------------------------------------------
# Obtaining the test data:

# Data source for 20-electrode raw runs:
RAW_E20_RUNS_SOURCE_DIR := ../../nmeeg_split_e20_gh2012/tests/out/

data/%_e20_r.txt: ${RAW_E20_RUNS_SOURCE_DIR}/%.txt
	cp -av ${RAW_E20_RUNS_SOURCE_DIR}/$*.txt data/$*_e20_r.txt

# Data source for 128-electrode raw runs:
RAW_E128_RUNS_SOURCE_DIR := ../../nmeeg_split_e128_gh2013/tests/out/

data/%_e128_r.txt: ${RAW_E128_RUNS_SOURCE_DIR}/%.txt
	cp -av ${RAW_E128_RUNS_SOURCE_DIR}/$*.txt data/$*_e128_r.txt

# Data source for filtered runs:
FILT_RUNS_SOURCE_DIR := ../../nmeeg_filter/tests/out/

data/%_f.txt: ${FILT_RUNS_SOURCE_DIR}/%_f.txt
	cp -av ${FILT_RUNS_SOURCE_DIR}/$*_f.txt data/$*_f.txt

# Script that builds synthetic runs:
SYNTH_EEG_MAKER := fake_eeg_data.gawk

data/%_s.txt: ${SYNTH_EEG_MAKER}
	./${SYNTH_EEG_MAKER} -v nt=3600 > data/$*_s.txt

# Script that builds impulse runs:
PULSE_EEG_MAKER := ${STOLFIHOME}/bin/nmeeg_make_test_pulses.gawk

data/%_i.txt: ${PULSE_EEG_MAKER}
	${PULSE_EEG_MAKER} -v varied=0 -v nt=3600 > data/$*_i.txt

# ======================================================================
endif
