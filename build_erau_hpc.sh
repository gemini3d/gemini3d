#!/bin/sh
#
# "-d" option makes this a Debug build
#
# This is for Matt's machine, so he can use the older libraries he's been working with for some time.

module load gcc/6.3.0
module load blas/gcc/64/3.7.0
module load lapack/gcc/64/3.7.0
#module load openmpi/gcc/64/1.10.3
module load openmpi/gcc/64/3.1.2
module load scalapack/openmpi/gcc/64/2.0.2

MUMPSPREFIX=$HOME/zettergmdata/lib/MUMPS_4.10.0

OPTS="-DNP=4"

export FC=/cm/shared/apps/openmpi/gcc/64/3.1.2/bin/mpifort
export CC=/cm/shared/apps/openmpi/gcc/64/3.1.2/bin/mpicc

. script_utils/check.sh

. script_utils/build.sh
