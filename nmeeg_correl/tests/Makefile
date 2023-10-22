# Last edited on 2021-08-26 08:59:46 by stolfi

PROG := nmeeg_correl
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# Format: {DATASET}:{SNAME}:{SVAL}:{NVAR}
# {SNAME} = name of channel to use as a frame selector, or "@" to use all frames.
# {SVAL} = value of marker that implies selection of a frame.
# {NVAR} = 1 to normalize variance.
 
TESTS := \
  S14_mfltB:@:0:0
  
E128_SETS := \
  S14_mfltB:@:0:0

SELECT_SEL_SETS := \
  s004_mrunA_y:SOn:30:0 \
  s004_mrunA_y:Occ:50:0 \
  s004_mrunA_y:RAp:70:0
  
MULTI_RUN_TESTS := \
  s004_mrunA_x:@:0:0 \
  s004_mrunA_y:@:0:0
  
LONG_TESTS := \
  s004_longA_x:@:0:0 \
  s004_longA_y:@:0:0 \
   \
  s006_longA_x:@:0:0 \
  s006_longA_y:@:0:0

# Process experimental runs
EXRUN_TESTS := \
  s001_r005_r:@:0:0 \
  s001_r039_r:@:0:0 \
   \
  s002_r016_r:@:0:0
  
FILTRUN_TESTS := \
  s001_r005_f:@:0:0 \
  s001_r039_f:@:0:0 \
   \
  s002_r016_f:@:0:0 \
  s002_r043_f:@:0:0 \
   \
  s003_r009_f:@:0:0 \
  s003_r030_f:@:0:0 \
   \
  s004_r030_f:@:0:0 \
  s004_r036_f:@:0:0 \
   \
  s005_r014_f:@:0:0 \
  s005_r044_f:@:0:0 \
   \
  s006_r024_f:@:0:0 \
  s006_r042_f:@:0:0 \
   \
  s007_r028_f:@:0:0 \
  s007_r046_f:@:0:0 \
   \
  s008_r015_f:@:0:0 \
  s008_r049_f:@:0:0

ALL_TESTS := \
  ${EXRUN_TESTS} \
  ${LONG_RAW_TESTS} \
  ${FILTRUN_TESTS}

.PHONY:: do-all-tests do-single-test do-all-shows do-single-show get-data

all: do-all-tests
# all: do-all-shows

do-all-tests:
	for tt in ${TESTS} ; do \
          ${MAKE} TEST=$$tt do-single-test do-single-show ; \
        done

do-all-shows:
	for tt in ${TESTS} ; do \
          ${MAKE} TEST=$$tt do-single-show ; \
        done

ifneq "/${TEST}" "/"
# ======================================================================
# Recursive section, client must define ${TEST}
  
TEST_ARGS := ${subst :, ,${TEST}}
DATASET :=  ${word 1,${TEST_ARGS}}
SNAME :=    ${word 2,${TEST_ARGS}}
SVAL :=     ${word 3,${TEST_ARGS}}
NVAR :=     ${word 4,${TEST_ARGS}}

SHOW := SHOW

DATA_DIR := data
OUT_DIR := out

EEG_FILE := ${DATA_DIR}/${DATASET}.txt

OUT_PREFIX := ${OUT_DIR}/${DATASET}_

ifeq "/${SNAME}" "/@"
  SELECT_OPT := 
else
  SELECT_OPT := -select ${SNAME} ${SVAL}
endif

do-single-test: ${PROGDIR}/${PROG} ${EEG_FILE}
	rm -f ${OUT_PREFIX}_*
	${PROGDIR}/${PROG} \
          -outPrefix ${OUT_PREFIX} \
          -normAvg 1 \
          -normVar ${NVAR} \
          ${SELECT_OPT} \
          -minMag 2.0 \
          < ${EEG_FILE}

do-single-show: ${OUT_PREFIX}P000_eig.txt
	if [[ "/${SHOW}" == "/SHOW" ]]; then \
          echo "plotting the electrode correlation as function of distance" 1>&2 ; \
          nmeeg_plot_dst_cor.sh SHOW ${OUT_PREFIX}dst_cov.txt ${OUT_PREFIX}dst_cov.png ; \
        fi
	echo "displaying the principal component patterns" 1>&2 ; \
          make_component_images.sh ${SHOW} ${OUT_PREFIX} ${OUT_PREFIX}P???_eig.txt
        
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
PROJECT_E20_DIR := projects/neuromat/00-DATA/2013-05-23-ghislain

${DATA_DIR}/s004_longA_x.txt: ${PROJECT_E20_DIR}/raw/s4_bl1.txt.gz ~/bin/nmeeg_make_big_data.sh
	nmeeg_make_big_data.sh 4 1 0 longA 30000 18000

${DATA_DIR}/s004_longA_y.txt: ${DATA_DIR}/s004_longA_x.txt ~/bin/nmeeg_make_big_data.sh
	nmeeg_make_big_data.sh 4 1 1 longA 30000 18000

${DATA_DIR}/s006_longA_x.txt: ${PROJECT_E20_DIR}/raw/s6_bl1.txt.gz ~/bin/nmeeg_make_big_data.sh
	nmeeg_make_big_data.sh 6 1 0 longA 10000 18000

${DATA_DIR}/s006_longA_y.txt: ${DATA_DIR}/s006_longA_x.txt ~/bin/nmeeg_make_big_data.sh
	nmeeg_make_big_data.sh 6 1 1 longA 10000 18000

# ----------------------------------------------------------------------
# Concatenation of many runs:
${DATA_DIR}/s004_mrunA_x.txt: ${PROJECT_E20_DIR}/raw-runs/s004/r001.txt ./make_mrun_data.sh
	./make_mrun_data.sh 004 0 10 40 mrunA

${DATA_DIR}/s004_mrunA_y.txt: ${PROJECT_E20_DIR}/flt-runs/s004/r001_f.txt ./make_mrun_data.sh
	./make_mrun_data.sh 004 1 10 40 mrunA

# ----------------------------------------------------------------------
# Concatenation of many runs (128 electrodes):
PROJECT_E128_DIR := projects/neuromat/00-DATA/2013-09-04-ghislain-128

${DATA_DIR}/S14_mfltB.txt: ${PROJECT_E128_DIR}/flt-runs-B/s014/r00101.txt 
	cat ${PROJECT_E128_DIR}/flt-runs-B/s014/r00{110,210,310,410,510,610,710,810}.txt  \
          > ${DATA_DIR}/S14_mfltB.txt

# ======================================================================
endif
	
