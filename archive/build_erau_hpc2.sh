#!/bin/bash
#
# "-d" option makes this a Debug build
# "-t" option makes this a Trace build (dump certain variables to disk)
#
# this convenience script initally (one-time) setups up Gemini for gfortran
# *** for subsequent builds, you can just type "make" in the objects/ directory ***
# (I keep a second Terminal tab for this purpose)

# This is for Matt's HPC

set -e
set -u

PREFIX=$HOME/zettergmdata/lib
SUFFIX=gcc7

module load gcc/6.3.0
module load blas/gcc/64/3.7.0
module load lapack/gcc/64/3.7.0
#module load openmpi/gcc/64/1.10.3
module load openmpi/gcc/64/3.1.2
module load scalapack/openmpi/gcc/64/2.0.2

OPTS="-DNP=4"

#======================================================
MPIPREFIX=/cm/shared/apps/openmpi/gcc/64/3.1.2/
LAPACKPREFIX=$PREFIX/lapack-$SUFFIX
SCALAPACKPREFIX=$PREFIX/scalapack-$SUFFIX
MUMPSPREFIX=$PREFIX/mumps-$SUFFIX

export FC=$MPIPREFIX/bin/mpifort
export CC=$MPIPREFIX/bin/mpicc

# ============================================================
# this temporarily disables Intel compiler (if installed) from messing up your gfortran environment.
MKLROOT=
LD_LIBRARY_PATH=

# ----- build ! -------
. script_utils/check.sh

. script_utils/build.sh

