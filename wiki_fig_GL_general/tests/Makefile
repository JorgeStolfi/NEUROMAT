# Last edited on 2019-10-19 03:05:06 by jstolfi


PROG := makefigs
PROGDIR = ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

OUTPREFIX := out/GL

all: do-test

do-test: do-clean ${PROGDIR}/${PROG}
	nice ${PROGDIR}/${PROG} -outName ${OUTPREFIX}
	atril ${OUTPREFIX}-*.eps

do-clean:
	rm -f ${OUTPREFIX}-*.eps
