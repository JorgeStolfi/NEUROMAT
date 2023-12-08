# Last edited on 2023-11-29 11:22:59 by stolfi

PROG := nmeeg_median_filter
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# Format: {DATASET}:{KEEP_MEAN}:{KIND_LO}:{PARM_LO}:{WD_LO}:{KIND_HI}:{PARM_HI}:{WD_HI}

# Some actual datasets
SOME_SYNTH_TESTS := \
  synth_p:1:bino:0.000:999:bino:0.000:099 \
  synth_p:1:gaus:099.0:751:gaus:009.9:075 \
  synth_p:1:hann:0.000:501:hann:0.000:051 \
   \
  synth_s:1:bino:0.000:999:bino:0.000:099 \
  synth_s:1:gaus:099.0:751:gaus:009.9:075 \
  synth_s:1:hann:0.000:501:hann:0.000:051
  
BAD_SYNTH_TESTS := \
  synth_p:1:hann:0.500:501:hann:0.500:051 \
  synth_p:1:hann:1.000:501:hann:1.000:051 \
   \
  synth_s:1:hann:0.500:501:hann:0.500:051 \
  synth_s:1:hann:1.000:501:hann:1.000:051 \
   
SYNTH_TESTS := \
  synth_i:1:hann:0.000:501:hann:0.000:051 \
  synth_m:1:hann:0.000:501:hann:0.000:051 \
  synth_p:1:hann:0.000:501:hann:0.000:051 \
  synth_s:1:hann:0.000:501:hann:0.000:051 \
  synth_s:1:hann:0.000:255:hann:0.000:051
 
SOME_REAL_TESTS := \
  s013_r101_e128_r:0:bino:0.000:999:bino:0.000:099 \
  s013_r101_e128_r:0:gaus:099.0:751:gaus:009.9:075 \
  s013_r101_e128_r:0:hann:0.000:501:hann:0.000:051
  
BAD_REAL_TESTS := \
  s013_r101_e128_r:0:hann:0.500:501:hann:0.500:051 \
  s013_r101_e128_r:0:hann:1.000:501:hann:1.000:051 \
  s013_r101_e128_r:0:tria:0.000:501:tria:0.000:051  

REAL_TESTS := \
  s001_r005_e020_r:0:hann:0.000:501:hann:0.000:051 \
  s003_r009_e020_r:0:hann:0.000:501:hann:0.000:051 \
  s005_r044_e020_r:0:hann:0.000:501:hann:0.000:051 \
  s013_r101_e128_r:0:hann:0.000:501:hann:0.000:051

# TESTS := ${SYNTH_TESTS} ${REAL_TESTS}
# TESTS := ${REAL_TESTS}

TESTS := ${SOME_SYNTH_TESTS}
# TESTS := ${SOME_REAL_TESTS}

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
DATASET :=   ${word 1,${TEST_ARGS}}
KEEP_MEAN := ${word 2,${TEST_ARGS}}
KIND_LO :=   ${word 3,${TEST_ARGS}}
PARM_LO :=   ${word 4,${TEST_ARGS}}
WD_LO :=     ${word 5,${TEST_ARGS}}
KIND_HI :=   ${word 6,${TEST_ARGS}}
PARM_HI :=   ${word 7,${TEST_ARGS}}
WD_HI :=     ${word 8,${TEST_ARGS}}

SHOW := SHOW

IPREF_EEG := data/${DATASET}

# Prefix for output unprocessed and processed files:
OPREF := out/${DATASET}_${KIND_LO}_${subst .,,${PARM_LO}}_${WD_LO}_${KIND_HI}_${subst .,,${PARM_HI}}_${WD_HI}

# Prefix for output files: original, filtered, noise, trend:
OPREF_MM := ${OPREF}_mm
OPREF_HI := ${OPREF}_hi
OPREF_LO := ${OPREF}_lo

IFILE_EEG := ${IPREF_EEG}.txt
OFILE_MM := ${OPREF_MM}.txt
OFILE_LO := ${OPREF_LO}.txt
OFILE_HI := ${OPREF_HI}.txt

NE := ${shell  head -n 200 ${IFILE_EEG} | grep 'ne = ' | gawk '//{ print $$3; }'}
CHANNELS := ${shell head -n 200 ${IFILE_EEG} | grep 'channels = ' | gawk -v ne=${NE} '//{ for (k=0; (k<ne)&&(k<5); k++){ printf " %s", $$(3+k); } }'}

# Max voltage (uV) for plot scales
VMAX := 160

do-single-test: ${PROGDIR}/${PROG} out ${IFILE_EEG}
	rm -f ${OUT_DAT_PREFIX}_* ${OPREF}_*
	@echo " "
	@echo "== median filtering of ${IPREF_EEG} widths ${WD_LO} ${WD_HI} =="
	@echo "NE = ${NE}"
	@echo "CHANNELS = ${CHANNELS}"
	@echo "== plotting the input file ${IFILE_EEG} =="
	nmeeg_plot_channels.sh   ${SHOW} ${IPREF_EEG} 0 ${VMAX}  0 9999  0 0  1400 500
	@echo " "
	${PROGDIR}/${PROG} \
          -lowWeights ${WD_LO} ${KIND_LO} ${PARM_LO} \
          -highWeights ${WD_HI} ${KIND_HI} ${PARM_HI} \
          -keepMean ${KEEP_MEAN} \
          -outPrefix ${OPREF} \
          -writeLow -writeHigh \
          -verbose \
          < ${IFILE_EEG}
	@echo "== plotting the filtered output ${OFILE_MM} =="
	nmeeg_plot_channels.sh   ${SHOW} ${OPREF_MM} 0 ${VMAX}  0 9999  0 0  1400 500
	@echo " "
	@echo "== plotting the low-freq componnet ${OFILE_LO} =="
	nmeeg_plot_channels.sh   ${SHOW} ${OPREF_LO} 0 ${VMAX}  0 9999  0 0  1400 500
	@echo " "
	@echo "== plotting the high-freq component ${OFILE_HI} =="
	nmeeg_plot_channels.sh   ${SHOW} ${OPREF_HI} 0 ${VMAX}  0 9999  0 0  1400 500
	@echo " "
        
# ======================================================================
endif

# ----------------------------------------------------------------------
# Obtaining the test data:

# ----------------------------------------------------------------------
# Data source for 20-electrode raw runs:
DATA_E20_SOURCE_DIR := projects/neuromat/00-DATA/2013-05-23-ghislain 
RAW_E20_RUNS_SOURCE_DIR := ${DATA_E20_SOURCE_DIR}/raw-runs

data/%_e20_r.txt: ${RAW_E20_RUNS_SOURCE_DIR}/%.txt
	cp -av ${RAW_E20_RUNS_SOURCE_DIR}/$*.txt data/$*_e20_r.txt

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
