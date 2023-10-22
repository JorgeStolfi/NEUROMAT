# Last edited on 2014-01-23 02:56:01 by stolfilocal

PROG := nmeeg_make_imgbasis
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# Format: {NX}:{NY}:{NE}:{BTYPE}
 
TESTS := \
  256:256:128:03
 
E20_TESTS := \
  280:320:20:00 \
  280:320:20:01 \
  280:320:20:02
 
E128_TESTS := \
  280:320:128:00 \
  280:320:128:01 \
  280:320:128:02

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
BTYPE := ${word 4,${TEST_ARGS}}

SHOW := NOSHOW

OUT_PREFIX := out/${NX}_${NY}_e${NE}

do-single-test: ${PROGDIR}/${PROG} out
	rm -f ${OUT_PREFIX}_*
	${PROGDIR}/${PROG} \
          -type ${BTYPE} \
          -size ${NX} ${NY} \
          -electrodes ${NE} \
          -outPrefix ${OUT_PREFIX}
	if [[ "${SHOW}" == "SHOW" ]]; then \
          for ff in ${OUT_PREFIX}_msk.fni ${OUT_PREFIX}_elp.fni ${OUT_PREFIX}_b${BTYPE}_*.fni ; do \
            fni_view -scale 50 $$ff ; \
          done ; \
        fi
        
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
