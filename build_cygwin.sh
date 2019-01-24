#!/bin/bash
#
# "-d" option makes this a Debug build
# "-t" option makes this a Trace build (dump certain variables to disk)
#
# this convenience script initally (one-time) setups up Gemini for gfortran
# *** for subsequent builds, you can just type "make" in the objects/ directory ***
# (I keep a second Terminal tab for this purpose)

set -e
set -u

PREFIX=/cygdrive/c/Users/ggrubbs1/Desktop/lib


#======================================================
MUMPSPREFIX=$PREFIX/MUMPS_5.1.1

export FC=mpifort
export CC=mpicc

# ============================================================
# this temporarily disables Intel compiler (if installed) from messing up your gfortran environment.
MKLROOT=
LD_LIBRARY_PATH=

for d in $MUMPSPREFIX
do
  [[ -d $d ]] || { echo "ERROR: $d not found"; exit 1; }
done

OPTS="-DMUMPS_ROOT=$MUMPSPREFIX ${OPTS:-}"

cmake --version

[[ ${1:-} == "-d" ]] && OPTS="-DCMAKE_BUILD_TYPE=Debug $OPTS"
[[ ${1:-} == "-t" ]] && OPTS="-DTRACE:BOOL=on $OPTS"

rm -rf objects/*  # need this one-time in case different compiler e.g. ifort was previously used.

cmake $OPTS -B objects -S .

cmake --build objects -j
