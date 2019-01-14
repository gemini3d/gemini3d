#!/bin/bash
# this convenience script initally (one-time) setups up Gemini for gfortran
# *** for subsequent builds, you can just type "make" in the objects/ directory ***
# (I keep a second Terminal tab for this purpose)

OPTS="-DMUMPS_ROOT=$HOME/flibs-gnu-mkl/MUMPS"

[[ $1 == "-d" ]] && OPTS="-DCMAKE_BUILD_TYPE=Debug $OPTS"
[[ $1 == "-t" ]] && OPTS="-DTRACE:BOOL=on $OPTS"

# if MKLROOT is not defined, try a default value
[[ -v MKLROOT ]] || export MKLROOT=$HOME/intel/compilers_and_libraries/linux/mkl/

. $MKLROOT/../bin/compilervars.sh intel64
. $MKLROOT/bin/mklvars.sh intel64 ilp64

export FC=/usr/bin/mpif90

rm -rf objects/*  # need this one-time in case different compiler e.g. ifort was previously used.

cmake $OPTS -B objects .

cmake --build objects -j

