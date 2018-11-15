#!/bin/bash
# this convenience script initally (one-time) setups up Gemini for gfortran
# *** for subsequent builds, you can just type "make" in the objects/ directory ***
# (I keep a second Terminal tab for this purpose)

# this temporarily disables Intel compiler (if installed) from messing up your gfortran environment.
MKLROOT=
FC=~/.local/bin/mpifort

# MUMPS is provided for Gfortran by:
# apt install libmumps-dev

#---exe---
(
rm -r objects/*  # need this one-time in case different compiler e.g. ifort was previously used.
cd objects
FC=~/.local/bin/mpifort cmake .. -DMUMPS_ROOT=/home/ggrubbs1/source/MUMPS_5.1.1 -DSCALAPACK_ROOT=/home/ggrubbs1/source/scalapack-2.0.2/ -DLAPACK_ROOT=/home/ggrubbs1/source/lapack-3.8.0/ -DBLAS_ROOT=/home/ggrubbs1/source/lapack-3.8.0/
)

cmake --build objects

