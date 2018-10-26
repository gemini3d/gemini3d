#!/bin/bash
# prereqs -- you'll need these in your current runtime shell as well to successfully run the program.
# *** for subsequent builds, you can just type "make" in the objects/ directory ***
# (I keep a second Terminal tab for this purpose)
# I.e. copy and paste them or put in ~/.bashrc
. $MKLROOT/../bin/compilervars.sh intel64
. $MKLROOT/bin/mklvars.sh intel64 ilp64
#------ MUMPS rebuild ----------------------------------------------------------------------
FC=mpiifort CC=mpiicc make d -s -C $MUMPS_ROOT -j2

#----- gemini -------------
(
rm -r objects/* # one-time, if you build for Gfortran previously
cd objects

FC=mpiifort cmake -DMUMPS_ROOT=$MUMPS_ROOT ..
)

cmake --build objects -j

