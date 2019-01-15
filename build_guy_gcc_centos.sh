#!/bin/bash
# this convenience script initally (one-time) setups up Gemini for gfortran
# *** for subsequent builds, you can just type "make" in the objects/ directory ***
# (I keep a second Terminal tab for this purpose)

set -e

cmake --version

# this temporarily disables Intel compiler (if installed) from messing up your gfortran environment.
export MKLROOT=
export FC=~/.local/bin/mpifort

# MUMPS is provided for Gfortran by:
# apt install libmumps-dev

rm -rf objects/*  # need this one-time in case different compiler e.g. ifort was previously used.

cmake -DMUMPS_ROOT=/home/ggrubbs1/source/MUMPS_5.1.1 -DSCALAPACK_ROOT=/home/ggrubbs1/source/scalapack-2.0.2/ -DLAPACK_ROOT=/home/ggrubbs1/source/lapack-3.8.0/ -DBLAS_ROOT=/home/ggrubbs1/source/lapack-3.8.0/ -B objects -S .

cmake --build objects -j

