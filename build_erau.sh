#!/bin/bash
#
# "-d" option makes this a Debug build
# "-t" option makes this a Trace build (dump certain variables to disk)
#
# This is for Matt's machine, so he can use the older libraries he's been working with for some time.


set -e
set -u

PREFIX=$HOME/lib
SUFFIX=gcc8

#======================================================
MPIPREFIX=/usr/lib64/openmpi/
LAPACKPREFIX=/usr/lib64/
SCALAPACKPREFIX=$PREFIX/scalapack-2.0.2
MUMPSPREFIX=$PREFIX/MUMPS_4.10.0

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

cmake --version

[[ ${1:-} == "-d" ]] && OPTS="-DCMAKE_BUILD_TYPE=Debug $OPTS"
[[ ${1:-} == "-t" ]] && OPTS="-DTRACE:BOOL=on $OPTS"

rm -rf objects/*  # need this one-time in case different compiler e.g. ifort was previously used.

cmake $OPTS -B objects -S .

cmake --build objects -j

