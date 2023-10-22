# Last edited on 2021-08-29 01:08:54 by stolfi

PROG := nmeeg_animate
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# Format: {DATASET}:{NE}:{N_SKIP}:{N_USE}:{RSTEP}:{B_TYPE}:{NX}:{NY}
 
TESTS := \
  s013_r00112_fC:129:000000:000700:1:voronoi_i0_n0:0450:0630

SINGLE_FRAME_TESTS := \
  s013_r00112_fC:129:000100:000001:1:voronoi_i0_n0:0450:0630 \
  s013_r00112_fC:129:000100:000001:1:shepard_i0_n1:0450:0630 \
  s013_r00112_fC:129:000100:000001:1:gauss_i1_n1:0450:0630

FILTERED_RUN_TESTS := \
  s013_r00112_fC:129:000050:000650:1:voronoi_i0_n0:0450:0630 \
  s013_r00112_fC:129:000050:000650:1:shepard_i0_n1:0450:0630 \
  s013_r00112_fC:129:000050:000650:1:gauss_i1_n1:0450:0630

RAW_RUN_TESTS := \
  s013_r00112_rw:20:000000:003800:5:voronoi_i0_n0:0450:0630 \
  s013_r00112_rw:20:001000:001800:1:voronoi_i0_n0:0450:0630 \

ALL_TESTS := \
  ${SINGLE_FRAME_TESTS} \
  ${FILTERED_RUN_TESTS} \
  ${RAW_RUN_TESTS} \

# ----------------------------------------------------------------------
# Test data sources:
#
#  s013_r00112_fC.txt = projects/neuromat/00-${IN_DIR}/2013-09-04-ghislain-128/flt-runs-C/s013/r00112.txt
#  s013_r00112_rr.txt = projects/neuromat/00-${IN_DIR}/2013-09-04-ghislain-128/raw-runs/s013/r00112.txt
# ----------------------------------------------------------------------

.PHONY:: \
  do-all-tests do-single-test \
  do-single-clean do-single-frames do-single-mp4 do-single-gif

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

SHOW := YES

IN_DIR := data

EEG_FILE := ${IN_DIR}/${DATASET}.txt

MSK_DIR := ${IN_DIR}/imgbasis/${NX}_${NY}_e${NE}
BAS_DIR := ${MSK_DIR}/${B_TYPE}

OUT_DIR := out/${NX}_${NY}_e${NE}/${B_TYPE}/${N_SKIP}_${N_USE}_s${RSTEP}
FRAME_DIR := ${OUT_DIR}/${DATASET}

GIF_FILE := ${FRAME_DIR}/movie.gif
MP4_FILE := ${FRAME_DIR}/movie.mp4

# do-single-test: do-single-clean do-single-frames do-single-gif
do-single-test: do-single-clean do-single-frames do-single-mp4
# do-single-test: do-single-mp4

do-single-clean:
	rm -f ${FRAME_DIR}/*.png ${GIF_FILE} ${MP4_FILE}

do-single-frames: ${PROGDIR}/${PROG} ${IN_DIR}/${DATASET}.txt
	mkdir -p ${FRAME_DIR}
	${PROGDIR}/${PROG} \
          -electrodes ${NE} \
          -basisDir ${BAS_DIR} \
          -maskDir ${MSK_DIR} \
          -skip ${N_SKIP} -read ${N_USE} -step ${RSTEP} \
          -inDir ${IN_DIR} \
          -outDir ${OUT_DIR} \
          -marker 129 1.00 1.000 0.150 0.000 \
          -marker 130 1.00 0.000 0.850 1.000 \
          -setName ${DATASET}

do-single-mp4: 
	do_ffmpeg.sh ${FRAME_DIR} ${MP4_FILE}

do-single-gif: 
	convert \
          ${MSK_DIR}/msk.png ${MSK_DIR}/msk.png ${MSK_DIR}/msk.png \
          `ls ${FRAME_DIR}/f??????.png | sort` -loop 2 -delay 5/100 \
          ${GIF_FILE}
	if [[ "${SHOW}" == "YES" ]]; then nomacs ${GIF_FILE}; fi
        

# ======================================================================
endif
