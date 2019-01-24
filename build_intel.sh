#!/bin/bash
#
# "-d" option makes this a Debug build
# "-t" option makes this a Trace build (dump certain variables to disk)
#
# this convenience script initally (one-time) setups up Gemini for ifort
# *** for subsequent builds, you can just type "make" in the objects/ directory ***
# (I keep a second Terminal tab for this purpose)

set -e
set -u

PREFIX=$HOME/.local
SUFFIX=ifort19

#======================================================
MPIPREFIX=$MKLROOT/../mpi/intel64
MUMPSPREFIX=$PREFIX/mumps-$SUFFIX


# some systems don't have mpiifort for Intel
export FC=$MPIPREFIX/bin/mpiifort
export CC=$MPIPREFIX/bin/mpiicc
export CXX=icpc
# ============================================================

for d in $MPIPREFIX $MUMPSPREFIX
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

