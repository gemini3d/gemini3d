#!/bin/bash
# this convenience script initally (one-time) setups up Gemini for gfortran
# *** for subsequent builds, you can just type "make" in the objects/ directory ***
# (I keep a second Terminal tab for this purpose)

set -e

cmake --version

# this temporarily disables Intel compiler (if installed) from messing up your gfortran environment.
export MKLROOT=
export FC=mpifort


rm -rf objects/*  # need this one-time in case different compiler e.g. ifort was previously used.

cmake -DMPI_Fortran_COMPILER=mpifort -DMUMPS_ROOT=/cygdrive/c/Users/ggrubbs1/Desktop/lib/MUMPS_5.1.1 -B objects -S .

cmake --build objects -j
