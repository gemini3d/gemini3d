#!/bin/sh
#
# "-d" option makes this a Debug build
#
# This is for Matt's machine, so he can use the older libraries he's been working with for some time.

set -e

cmake --version

module load gcc/6.3.0
module load blas/gcc/64/3.7.0
module load lapack/gcc/64/3.7.0
#module load openmpi/gcc/64/1.10.3
module load openmpi/gcc/64/3.1.2
module load scalapack/openmpi/gcc/64/2.0.2

OPTS="-DMUMPS_ROOT=$HOME/zettergmdata/lib/MUMPS_4.10.0 $OPTS"
OPTS="-DNP=4 $OPTS"

export FC=/cm/shared/apps/openmpi/gcc/64/3.1.2/bin/mpifort
export CC=/cm/shared/apps/openmpi/gcc/64/3.1.2/bin/mpicc

[[ $1 == "-d" ]] && OPTS="-DCMAKE_BUILD_TYPE=Debug $OPTS"
[[ $1 == "-t" ]] && OPTS="-DTRACE:BOOL=on $OPTS"

rm -rf objects/*  # need this one-time in case different compiler e.g. ifort was previously used.

cmake $OPTS -B objects -S .

cmake --build objects -j
