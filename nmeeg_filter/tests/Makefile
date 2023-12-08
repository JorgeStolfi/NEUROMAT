# Last edited on 2023-12-07 17:53:40 by stolfi

PROG := nmeeg_filter
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make
 
# Define ${TESTS} with "=" not ":="
# TESTS = ${ALL_TESTS}
# TESTS = ${PULSE_TESTS}
# TESTS = ${IMPULSE_TESTS} ${PULSE_TESTS} ${RAMP_TESTS}
# TESTS = ${RAMP_TESTS}
# TESTS = ${SINE_TESTS}
TESTS = ${E20_TESTS} ${E128_TESTS}
# TESTS = synth_m:G:3
# TESTS = synth_p:BUR8:0
# TESTS = synth_p:G:0

# Format: {DATASET}:{FILTER_TYPE}:{DEG_TREND}:{EXCLUDE_CHANS}

# 128-electrode runs:
E128_TESTS := \
  s013_r00101_e128_r:G:1:0 \
  s013_r00101_e128_r:G:2:0 \
  s013_r00101_e128_r:G:3:0 \
  \
  s013_r00202_e128_r:G:3:0 \
  s013_r00402_e128_r:G:3:0 \
  s013_r00602_e128_r:G:3:0

# 20-electrode runs:
E20_TESTS := \
  s003_r009_e20_r:G:3:0 \
  s003_r030_e20_r:G:3:0 \
   \
  s001_r005_e20_r:G:3:0 \
  s001_r039_e20_r:G:3:0 \
   \
  s002_r016_e20_r:G:3:0 \
  s002_r043_e20_r:G:3:0 \
   \
  s004_r030_e20_r:G:3:0 \
  s004_r036_e20_r:G:3:0 \
   \
  s005_r014_e20_r:G:3:0 \
  s005_r044_e20_r:G:3:0 \
   \
  s006_r024_e20_r:G:3:0 \
  s006_r042_e20_r:G:3:0 \
   \
  s007_r028_e20_r:G:3:0 \
  s007_r046_e20_r:G:3:0 \
   \
  s008_r015_e20_r:G:3:0 \
  s008_r049_e20_r:G:3:0

IMPULSE_TESTS := \
  synth_i:G:0:1 \
  synth_i:G:3:1
  
PULSE_TESTS := \
  synth_i:G:0:1 \
  synth_i:G:3:1

RAMP_TESTS := \
  synth_m:G:0:1 \
  synth_m:G:1:1 \
  synth_m:G:2:1 \
  synth_m:G:3:1
  
SINE_TESTS := \
  synth_s:G:0:1 \
  synth_s:G:1:1 \
  synth_s:G:2:1 \
  synth_s:G:3:1
  
FTYPE_TESTS :=  \
  synth_i:G:0 \
  synth_p:G:0 \
  \
  synth_i:BUR8:0 \
  synth_p:BUR8:0 \
  \
  synth_i:BUC12:0 \
  synth_p:BUC12:0 \
  \
  synth_i:SG:0 \
  synth_p:SG:0 \
  \
  synth_i:DL:0 \
  synth_p:DL:0 \
  \
  synth_i:BQ:0 \
  synth_p:BQ:0

ALL_TESTS := \
  ${PULSE_TESTS} \
  ${IMPULSE_TESTS} \
  ${RAMP_TESTS} \
  ${SINE_TESTS} \
  ${FTYPE_TESTS} \
  ${E20_TESTS} \
  ${E128_TESTS}

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
DATASET :=       ${word 1,${TEST_ARGS}}
FILTER_TYPE :=   ${word 2,${TEST_ARGS}}
DEG_TREND :=     ${word 3,${TEST_ARGS}}
EXCLUDE_CHANS := ${word 4,${TEST_ARGS}}

SHOW := SHOW

EEG_FILE := data/${DATASET}.txt

# Prefix for output unprocessed and processed files:
OUT_DAT_PREFIX := out/${DATASET}
OUT_PRC_PREFIX := out/${DATASET}_fty${FILTER_TYPE}_dgt${DEG_TREND}

# Prefix for output files: original, filtered, noise, trend:
ORG_PREFIX := ${OUT_DAT_PREFIX}
FLT_PREFIX := ${OUT_PRC_PREFIX}_ot
NSE_PREFIX := ${OUT_PRC_PREFIX}_ns
TRE_PREFIX := ${OUT_PRC_PREFIX}_tr

ORG_FILE := ${ORG_PREFIX}.txt
FLT_FILE := ${FLT_PREFIX}.txt
NSE_FILE := ${NSE_PREFIX}.txt
TRE_FILE := ${TRE_PREFIX}.txt

NE := ${shell cat ${EEG_FILE} | grep 'ne = ' | gawk '//{ print $$3; }'}
ifeq "/${EXCLUDE_CHANS}" "/1"
  # Test the exclude option:
  EXCLUDE_OPS := \
    -exclude C002 \
    -exclude ${LEX} \
    -exclude T01
else
  EXCLUDE_OPS := 
endif 
LEX := ${shell printf "C%03d" "${NE}"}

# Max voltage (uV) for plot scales
VMAX := 160

# Filter parameters, output resampling step, max freq to plot in resampled output:

# FLT_PARMS :=  ${FILTER_TYPE} 0.166667 0.333333  15.0 30.0
FLT_PARMS :=  ${FILTER_TYPE} 0.166667 0.333333  24.0 48.0
RSTEP := 5
FRSMAX := 300

do-single-test: ${PROGDIR}/${PROG} out ${EEG_FILE}
	rm -f ${OUT_DAT_PREFIX}_* ${OUT_PRC_PREFIX}_*
	@echo " "
	@echo "== plotting the input data and spectrum =="
	@echo "NE = ${NE}"
	nmeeg_plot_channels.sh   ${SHOW} data/${DATASET} 0 ${VMAX}  0 9999  0 0  1400 500
	nmeeg_spectrum           ${ORG_PREFIX} 0 0 < data/${DATASET}.txt
	nmeeg_plot_spectra.sh    ${SHOW} ${ORG_PREFIX}_pwr 300 0 $$(( ${NE} - 1 )) 480 'linespoints pt 7'
	@echo " "
	@echo "== filtering with trend preservation =="
	${PROGDIR}/${PROG} \
          ${EXCLUDE_OPS} \
          -filter ${FLT_PARMS} \
          -invert 0 \
          -resample ${RSTEP} \
          -trend ${DEG_TREND} 1 \
          -keepMean T \
          -verbose \
          < ${EEG_FILE} \
          > ${FLT_FILE}
	@echo "== plotting the filtered output and spectrum =="
	nmeeg_plot_channels.sh   ${SHOW} ${FLT_PREFIX} 0 ${VMAX}  0 9999  0 0  1400 500
	nmeeg_spectrum           ${FLT_PREFIX} 0 0 < ${FLT_FILE}
	nmeeg_plot_spectra.sh    ${SHOW} ${FLT_PREFIX}_pwr ${FRSMAX} 0 $$(( ${NE} - 3 )) 480 'linespoints pt 7'
	@echo " "
	@echo "== noise without trend =="
	${PROGDIR}/${PROG} \
          ${EXCLUDE_OPS} \
          -filter ${FLT_PARMS} \
          -invert 1 \
          -trend ${DEG_TREND} 0 \
          -keepMean T \
          -verbose \
          < ${EEG_FILE} \
          > ${NSE_FILE}
	@echo "== plotting the removed noise and spectrum =="
	nmeeg_plot_channels.sh   ${SHOW} ${NSE_PREFIX} 0 ${VMAX}  0 9999  0 0  1400 500
	nmeeg_spectrum           ${NSE_PREFIX} 0 0 < ${NSE_FILE}
	nmeeg_plot_spectra.sh    ${SHOW} ${NSE_PREFIX}_pwr 300 0 $$(( ${NE} - 1 )) 480 'linespoints pt 7'
	@echo " "
	@echo "== trend only =="
	${PROGDIR}/${PROG} \
          ${EXCLUDE_OPS} \
          -invert 1 \
          -trend ${DEG_TREND} 1 \
          -keepMean T \
          -verbose \
          < ${EEG_FILE} \
          > ${TRE_FILE}
	@echo "== plotting the trend only =="
	nmeeg_plot_channels.sh   ${SHOW} ${TRE_PREFIX} 0 ${VMAX}  0 9999  0 0  1400 500
	nmeeg_spectrum           ${TRE_PREFIX} 0 0 < ${TRE_FILE}
	nmeeg_plot_spectra.sh    ${SHOW} ${TRE_PREFIX}_pwr 300 0 $$(( ${NE} - 1 )) 480 'linespoints pt 7'
        
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
