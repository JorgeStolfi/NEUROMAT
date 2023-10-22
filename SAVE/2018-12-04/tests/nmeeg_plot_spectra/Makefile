# Last edited on 2013-06-13 23:49:59 by stolfilocal

PROG := nmeeg_plot_spectra.sh
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# Format: {DATASET}:{FMAX}:{STYLE}

TESTS := \
  s8_r015_f_pwr:080:popsicles \
  s8_r015_f_pwr:080:lines \
  s8_r015_f_pwr:080:linespoints \
  s8_r015_f_pwr:080:histeps \
   \
  s8_r015_r_pwr:300:lines

ALL_TESTS := \
  s8_r015_f_pwr:080:popsicles \
  s8_r015_f_pwr:080:lines \
  s8_r015_f_pwr:080:linespoints \
  s8_r015_f_pwr:080:histeps \
   \
  s8_r015_r_pwr:300:lines

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
FMAX :=     ${word 2,${TEST_ARGS}}
STYLE :=    ${subst _, ,${word 3,${TEST_ARGS}}}


SHOW := SHOW

EEG_FILE := data/${DATASET}.txt

OUT_PREFIX := out/${DATASET}_${N_SKIP}:${N_USE}

do-single-test: ${PROGDIR}/${PROG}
	rm -f ${OUT_PREFIX}-*
	${PROGDIR}/${PROG} \
          ${SHOW} \
          data/${DATASET} \
          ${FMAX} 0 19 500 '${STYLE}'

# ----------------------------------------------------------------------
# Obtaining the test data:

# Data source for raw runs:
RAW_RUNS_SOURCE_DIR := ../../nmeeg_split_e20_gh2012/tests/out/

data/%_r.txt: ${RAW_RUNS_SOURCE_DIR}/%.txt
	cp -av ${RAW_RUNS_SOURCE_DIR}/$*.txt data/$*_r.txt

# Data source for filtered runs:
FILT_RUNS_SOURCE_DIR := ../../nmeeg_filter/tests/out/

data/%_f.txt: ${FILT_RUNS_SOURCE_DIR}/%_r_000000_004800_f.txt
	cp -av ${FILT_RUNS_SOURCE_DIR}/$*_f.txt data/$*_f.txt

# ======================================================================
endif
