# Last edited on 2023-12-06 07:07:15 by stolfi

PROG := nmeeg_get_stats
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make
 
# Format: {DATASET}:{HSTEP}
TESTS := \
  synth_s:1.000 \
  s003_r009_e20_r:1.000

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
HSTEP :=  ${word 2,${TEST_ARGS}}

SHOW := SHOW

EEG_FILE := data/${DATASET}.txt

# Prefix for output histogram files:
OUT_PREFIX := out/${DATASET}_
OUT_STATS_FILE := ${OUT_PREFIX}stats.txt

NE := ${shell  head -n 200 ${EEG_FILE} | grep 'ne = ' | gawk '//{ print $$3; }'}
CHANNELS := ${shell head -n 200 ${EEG_FILE} | grep 'channels = ' | gawk -v ne=${NE} '//{ for (k=0; (k<ne)&&(k<5); k++){ printf " %s", $$(3+k); } }'}

do-single-test: ${PROGDIR}/${PROG} out ${EEG_FILE}
	rm -f ${OUT_PREFIX}*
	@echo " "
	@echo "== stats and histogram of ${DATASET} =="
	@echo "NE = ${NE}"
	${PROGDIR}/${PROG} \
          -verbose \
          -hStep ${HSTEP} \
          -outPrefix ${OUT_PREFIX} \
          < ${EEG_FILE} \
          > ${OUT_STATS_FILE}
	@echo "== plotting the histograms =="
	for ch in ${CHANNELS} ; do \
          hfile="${OUT_PREFIX}$${ch}_hist" ; \
          pfile="${OUT_PREFIX}$${ch}_hist" ; \
          nmeeg_plot_electrode_histogram.sh  ${SHOW} $${ch} $${hfile} $${pfile} ; \
        done
        
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# Obtaining the test data:

# ----------------------------------------------------------------------
# Data source for 20-electrode raw runs:
DATA_E20_SOURCE_DIR := projects/neuromat/00-DATA/2013-05-23-ghislain 
RAW_E20_RUNS_SOURCE_DIR := ${DATA_E20_SOURCE_DIR}/raw-runs

data/%_e20_r.txt: ${RAW_E20_RUNS_SOURCE_DIR}/%.txt
	cp -av ${RAW_E20_RUNS_SOURCE_DIR}/$*.txt data/$*_e20_r.txt

# ----------------------------------------------------------------------
# Raw data for 128-electrode raw runs:
DATA_E128_SOURCE_DIR := projects/neuromat/00-DATA/2013-09-04-ghislain-128
RAW_E128_RUNS_SOURCE_DIR := ${DATA_E128_SOURCE_DIR}/raw-runs

data/%_e128_r.txt: ${RAW_E128_RUNS_SOURCE_DIR}/%.txt
	cp -av ${RAW_E128_RUNS_SOURCE_DIR}/$*.txt data/$*_e128_r.txt

# ----------------------------------------------------------------------
# Synthetic EEG signals:
SYNTH_EEG_MAKER := ${STOLFIHOME}/bin/nmeeg_make_test_pulses.gawk

# Impulses:
data/synth_i.txt: ${SYNTH_EEG_MAKER}
	${SYNTH_EEG_MAKER} -v type=0 -v nt=3500 -v fsmp=500 > data/synth_i.txt

# Misc pulses:
data/synth_p.txt: ${SYNTH_EEG_MAKER}
	${SYNTH_EEG_MAKER} -v type=1 -v nt=3500 -v fsmp=500 > data/synth_p.txt

# Ramps:
data/synth_m.txt: ${SYNTH_EEG_MAKER}
	${SYNTH_EEG_MAKER} -v type=2 -v nt=3500 -v fsmp=500 > data/synth_m.txt

# Sinusoids:

data/synth_s.txt: ${SYNTH_EEG_MAKER}
	${SYNTH_EEG_MAKER} -v type=3 -v nt=3500 -v fsmp=500 > data/synth_s.txt

# ----------------------------------------------------------------------
# Raw data spanning multiple runs:
PROJECT_DIR := projects/neuromat/00-DATA/2013-05-23-ghislain

data/s004_longA_x.txt: ${PROJECT_DIR}/raw/s4_bl1.txt.gz ~/bin/nmeeg_make_big_data.sh
	nmeeg_make_big_data.sh 4 1 0 longA 30000 18000

data/s004_longA_y.txt: data/s004_longA_x.txt ~/bin/nmeeg_make_big_data.sh
	nmeeg_make_big_data.sh 4 1 1 longA 30000 18000

# ======================================================================
endif
