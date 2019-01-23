#!/bin/bash
#
# "-d" option makes this a Debug build
# "-t" option makes this a Trace build (dump certain variables to disk)
#
# this convenience script initally (one-time) setups up Gemini for gfortran
# *** for subsequent builds, you can just type "make" in the objects/ directory ***
# (I keep a second Terminal tab for this purpose)

set -e

cmake --version


OPTS="-DMPI_ROOT=$HOME/zettergmdata/lib-gcc7.3/openmpi-3.1.3 $OPTS"
OPTS="-DLAPACK_ROOT=$HOME/zettergmdata/lib-gcc7.3/lapack $OPTS"
OPTS="-DSCALAPACK_ROOT=$HOME/zettergmdata/lib-gcc7.3/fortran-libs/scalapack $OPTS"
OPTS="-DMUMPS_ROOT=$HOME/zettergmdata/lib-gcc7.3/fortran-libs/MUMPS $OPTS"
OPTS="-DNP=10 $OPTS"

export FC=$HOME/zettergmdata/lib-gcc7.3/openmpi-3.1.3/bin/mpifort
export CC=$HOME/zettergmdata/lib-gcc7.3/openmpi-3.1.3/bin/mpicc

# this temporarily disables Intel compiler (if installed) from messing up your gfortran environment.
MKLROOT=
LD_LIBRARY_PATH=

[[ $1 == "-d" ]] && OPTS="-DCMAKE_BUILD_TYPE=Debug $OPTS"
[[ $1 == "-t" ]] && OPTS="-DTRACE:BOOL=on $OPTS"

rm -rf objects/*  # need this one-time in case different compiler e.g. ifort was previously used.

cmake $OPTS -B objects -S .

cmake --build objects -j

