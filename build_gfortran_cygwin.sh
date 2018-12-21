#!/bin/bash
# this convenience script initally (one-time) setups up Gemini for gfortran
# *** for subsequent builds, you can just type "make" in the objects/ directory ***
# (I keep a second Terminal tab for this purpose)

# this temporarily disables Intel compiler (if installed) from messing up your gfortran environment.
MKLROOT=
FC=mpifort

# MUMPS is provided for Gfortran by:
# apt install libmumps-dev

#---exe---
(
rm -r objects/*  # need this one-time in case different compiler e.g. ifort was previously used.
cd objects
FC=mpifort cmake .. -DMPI_Fortran_COMPILER=mpifort -DMUMPS_ROOT=/cygdrive/c/Users/ggrubbs1/Desktop/lib/MUMPS_5.1.1
)

cmake --build objects
