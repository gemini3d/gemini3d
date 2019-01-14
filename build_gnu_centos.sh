#!/bin/bash
#
# "-d" option makes this a Debug build
# "-t" option makes this a Trace build (dump certain variables to disk)
#
# this convenience script initally (one-time) setups up Gemini for gfortran
# *** for subsequent builds, you can just type "make" in the objects/ directory ***
# (I keep a second Terminal tab for this purpose)

#OPTS="-DMPI_ROOT=$HOME/.local/openmpi-3.1-gcc7 $OPTS"
#OPTS="-DLAPACK_ROOT=$HOME/.local/lapack-gcc7 $OPTS"
OPTS="-DSCALAPACK_ROOT=$HOME/flibs-gnu-nomkl/scalapack $OPTS"
OPTS="-DMUMPS_ROOT=$HOME/flibs-gnu-nomkl/MUMPS $OPTS"

# this temporarily disables Intel compiler (if installed) from messing up your gfortran environment.
MKLROOT=
LD_LIBRARY_PATH=

[[ $1 == "-d" ]] && OPTS="-DCMAKE_BUILD_TYPE=Debug $OPTS"
[[ $1 == "-t" ]] && OPTS="-DTRACE:BOOL=on $OPTS"

export FC=mpif90

rm -rf objects/*  # need this one-time in case different compiler e.g. ifort was previously used.

cmake $OPTS -B objects .

cmake --build objects -j

