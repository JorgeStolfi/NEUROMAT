# Last edited on 2023-12-08 08:39:52 by stolfi


PROG := makefigs
PROGDIR = ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

OUTPREFIX := out/GL

all: do-test

do-test: do-clean ${PROGDIR}/${PROG}
	nice ${PROGDIR}/${PROG} -scale 10.0 -outName ${OUTPREFIX}
	atril ${OUTPREFIX}_*.eps

do-clean:
	rm -f ${OUTPREFIX}-*.eps
