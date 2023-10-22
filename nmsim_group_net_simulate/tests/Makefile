# Last edited on 2019-06-11 11:38:31 by jstolfi

PROG := nmsim_group_net_simulate
PROGDIR = ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

all: do-test

NETNAME := single
PARMSTAG := 01

TMIN := 50
TMAX := 550

NETFILE := in/${NETNAME}.txt
PARMSFILE := in/${NETNAME}_p${PARMSTAG}.parms

OUTPREFIX := out/${NETNAME}_p${PARMSTAG}

do-test: do-clean ${PROGDIR}/${PROG}
	nice ${PROGDIR}/${PROG} \
	  -netFile ${NETFILE} \
	  -timeStep 1.0 \
	  ` cat ${PARMSFILE} ` \
	  -outPrefix ${OUTPREFIX}
	for ff in ${OUTPREFIX}_*.txt ; do \
	  ./plot_potentials.sh ${TMIN} ${TMAX} $$ff ; \
	  ./plot_input_tot.sh ${TMIN} ${TMAX} $$ff ; \
	  ./plot_input_ext.sh ${TMIN} ${TMAX} $$ff ; \
	done

do-clean:
	rm -f ${OUTPREFIX}_*.txt ${OUTPREFIX}_*.png
