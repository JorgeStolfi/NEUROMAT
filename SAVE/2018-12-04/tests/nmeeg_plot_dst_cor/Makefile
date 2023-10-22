# Last edited on 2013-10-02 04:53:06 by stolfilocal

PROG := nmeeg_plot_dst_cor.sh
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

DATA_DIR := data
OUT_DIR := out

DATASET := test
FILE_NAME := ${DATASET}_dst_cov

SHOW := SHOW
OUT_PREFIX := ${OUT_DIR}/${DATASET}_dst_cov

DST_COV_FILE := ${DATA_DIR}/${FILE_NAME}.txt
PLOT_FILE := ${OUT_DIR}/${FILE_NAME}.png

all: ${PROGDIR}/${PROG} ${DST_COV_FILE}
	rm -f ${OUT_PREFIX}_*
	${PROGDIR}/${PROG} \
          ${SHOW} \
          ${DATA_DIR} ${OUT_DIR} ${FILE_NAME}
