# Last edited on 2023-10-22 19:41:07 by stolfi

IGNOREDIRS := \
  nmeeg_activity_plot \
  nmsim_mu0_diagram \
  nmsim_multipop
  
include GENERIC-ROOT-DIR.make

.PHONY:: find-scripts

all: uninstall build-progs install 

find-scripts:
	find ./ ~/projects/neuromat \( -name Makefile -o -name '*.sh' -o -name '00*.txt' \) > .scripts
