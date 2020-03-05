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

builddir=build

PREFIX=$HOME/zettergmdata/lib-gcc7.3

MPIPREFIX=$PREFIX/openmpi-3.1.3
LAPACKPREFIX=$PREFIX/lapack
SCALAPACKPREFIX=$PREFIX/fortran-libs/scalapack
MUMPSPREFIX=$PREFIX/fortran-libs/MUMPS

export FC=$MPIPREFIX/bin/mpifort
export CC=$MPIPREFIX/bin/mpicc


# ============================================================
# this temporarily disables Intel compiler (if installed) from messing up your gfortran environment.
MKLROOT=
LD_LIBRARY_PATH=

for d in $MPIPREFIX $LAPACKPREFIX $SCALAPACKPREFIX $MUMPSPREFIX
do
  [[ -d $d ]] || { echo "ERROR: $d not found"; exit 1; }
done

OPTS="-DMPI_ROOT=$MPIPREFIX ${OPTS:-}"
OPTS="-DLAPACK_ROOT=$LAPACKPREFIX $OPTS"
OPTS="-DSCALAPACK_ROOT=$SCALAPACKPREFIX $OPTS"
OPTS="-DMUMPS_ROOT=$MUMPSPREFIX $OPTS"
OPTS="-DNP=10 $OPTS"

cmake --version

[[ ${1:-} == "-d" ]] && OPTS="-DCMAKE_BUILD_TYPE=Debug $OPTS"
[[ ${1:-} == "-t" ]] && OPTS="-DTRACE:BOOL=on $OPTS"

cmake -B $builddir $OPTS

cmake --build $builddir -j

