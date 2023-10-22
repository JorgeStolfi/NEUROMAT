# nmsim_elem_net_simulate

## Attributes:

  + **Authors and collaborators:** Jorge Stolfi (IC-UNICAMP)
  + **Supervision:** Jorge Stolfi
  + **Intended users:** Neuroscience researchers
  + **Created:** Jorge Stolfi (2018)
  + **Situation:** Working, under development
  + **Info updated on:** 2014-04
  
## Dexcription and purpose:

This package contains a C program `nmsim_elem_net_simulate.c` that simulates a neuronal network using the Galves-Loecherbach neuron model, under the synchronous coarse discrete-time approach (time step ~1ms, no multi-step synapse or spike profile).

This program is similar to, but independent from, the synchronous discrete-time simulator developed by Nilton Kamiji at USP-RP, and from the simulator used by Antonio Roque also at USP-RP, which can also do fine-scale simulations (e.g. with time step 0.1 ms).

This program also differs from the continuos-time simulator being developed by Aline D. Oliveira at USP-IME-NUMEC.

  Folders:
  
    + **tests** - scripts and `Makefile` to run tests.
    
    + **tests/in** - some test networks and parameters.
    
    + **tests/out** - outputs of the test runs.
    
## Installation:

The top-level `Makefile` should compile the simulator `nmsim_elem_net_simulate.c`, but must be edited to fit your local conditions. 

The program requires the libraries `libnmsim_e.a`, `libnmsim.a`, and `jslibs.a` and the corresponding C header files, which can be obtained from the GitHub package `https://github.com/JorgeStolfi/JSLIBS`.
  
## Use

Execute `nmsim_elem_net_simulate --info` to get the full manpage of the program, and `nmsim_elem_net_simulate --help` for a short synopsis of the command line options.  

Basically the program reads a detailed network description from a ".txt" file, and simulates it over a specified time interval.  The external inputs and trace information are specified on the command line.  See the file `tests/Makefile` for examples.


Last edited on 2021-04-17 23:41:51 by jstolfi
