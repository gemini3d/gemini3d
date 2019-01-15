#!/bin/bash
#
# "-d" option makes this a Debug build
# "-t" option makes this a Trace build (dump certain variables to disk)
#
# *** for subsequent builds, you can just type "make" in the objects/ directory ***
# (I keep a second Terminal tab for this purpose)
# I.e. copy and paste them or put in ~/.bashrc

# --- if you haven't already setup your environment, you will need to do so for future "make" after this script

# if MKLROOT is not defined, try a default value
set -e

cmake --version

[[ -z ${MKLROOT+x} ]] && { echo 'MKLROOT must be set to use Intel compilers'; exit 1; }

# bash >= 4.2, centos 7 has bash 4.1
#[[ -v MKLROOT ]] || export MKLROOT=$HOME/intel/compilers_and_libraries/linux/mkl/

. $MKLROOT/../bin/compilervars.sh intel64
. $MKLROOT/bin/mklvars.sh intel64 ilp64

# some HPC installations have non-standard paths and names. Watch the beginning of the CMake output to be sure your Intel compiler is picked.
export FC=mpifort
export CC=mpicc
export CXX=icpc

rm -rf objects/* # one-time, if you build for Gfortran previously


# some systems don't have mpiifort for Intel
# use ifort as mpif90 get partially picked-up as GNU
OPTS="-DMUMPS_ROOT=$HOME/flibs-intel/MUMPS"

[[ $1 == "-d" ]] && OPTS="-DCMAKE_BUILD_TYPE=Debug $OPTS"
[[ $1 == "-t" ]] && OPTS="-DTRACE:BOOL=on $OPTS"

cmake $OPTS -B objects -S .

cmake --build objects -j

