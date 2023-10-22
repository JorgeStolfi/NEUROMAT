# Last edited on 2021-08-28 02:48:18 by stolfi

PROG := nmeeg_make_imgbasis
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# Format: {NX}:{NY}:{NE}:{BTYPE}:{INTERP}:{NORM}:{MSUB}

TNX := 0450
TNY := 0630

TNE = 129
 
TESTS = \
  ${TNX}:${TNY}:${TNE}:shepard:0:1:1
 
ALL_TESTS := \
  ${TNX}:${TNY}:${TNE}:shepard:0:1:1 \
  ${TNX}:${TNY}:${TNE}:gauss:0:0:1 \
  ${TNX}:${TNY}:${TNE}:gauss:0:1:1 \
  ${TNX}:${TNY}:${TNE}:gauss:1:0:1 \
  ${TNX}:${TNY}:${TNE}:gauss:1:1:1 \
  ${TNX}:${TNY}:${TNE}:mexican:0:0:1 \
  ${TNX}:${TNY}:${TNE}:mexican:0:1:1 \
  ${TNX}:${TNY}:${TNE}:mexican:1:1:1 \
  ${TNX}:${TNY}:${TNE}:mexican:1:0:1 \
  ${TNX}:${TNY}:${TNE}:voronoi:0:0:2

SHOW_PNG := YES
SHOW_FNI := NO

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
NX :=         ${word 1,${TEST_ARGS}}
NY :=         ${word 2,${TEST_ARGS}}
NE :=         ${word 3,${TEST_ARGS}}
BTYPE :=      ${word 4,${TEST_ARGS}}
INTERP :=     ${word 5,${TEST_ARGS}}
NORM :=       ${word 6,${TEST_ARGS}}
MSUB :=       ${word 7,${TEST_ARGS}}

OUTDIR := out
IMGDIR := ${OUTDIR}/${NX}_${NY}_e${NE}
BASDIR := ${IMGDIR}/${BTYPE}_i${INTERP}_n${NORM}

do-single-test: ${PROGDIR}/${PROG} out
	mkdir -p ${BASDIR}
	rm -f {${OUTDIR},${BASDIR}}/*.{fni,png}
	${PROGDIR}/${PROG} \
          -type ${BTYPE} \
          -interp ${INTERP} \
          -norm ${NORM} \
          -scalpSize 150 180 150 \
          -imageSize ${NX} ${NY} \
          -extent 0.95 \
          -electrodes ${NE} \
          -subsample ${MSUB} \
          -outDir ${OUTDIR}
	for fname in \
          ${IMGDIR}/msk.fni \
          ${IMGDIR}/elp.fni \
          ${BASDIR}/*.fni \
        ; do \
          pname="$${fname%.*}.png" ; \
          fni_to_png \
              -isMask T \
              -min -2.0 -max +2.0 \
              -maxval 65535 \
              -nanval 0 \
            < $${fname} \
            > $${pname} ; \
          if [[ "${SHOW_FNI}" == "YES" ]]; then fni_view -scale 50 $${fname} ; fi ; \
        done
	if [[ "${SHOW_PNG}" == "YES" ]]; then \
          eom ${OUTDIR}/msk.png ${OUTDIR}/elp.png ${BASDIR}/*.png; \
        fi ; \

out:
	mkdir out

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
