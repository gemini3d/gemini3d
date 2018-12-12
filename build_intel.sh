#!/bin/bash
#
# "-d" option makes this a Debug build
#
# *** for subsequent builds, you can just type "make" in the objects/ directory ***
# (I keep a second Terminal tab for this purpose)
# I.e. copy and paste them or put in ~/.bashrc

# --- if you haven't already setup your environment, you will need to do so for future "make" after this script

# if MKLROOT is not defined, try a default value
[[ -v MKLROOT ]] || export MKLROOT=$HOME/intel/compilers_and_libraries/linux/mkl/

. $MKLROOT/../bin/compilervars.sh intel64
. $MKLROOT/bin/mklvars.sh intel64 ilp64

# DO NOT change to mpif90 or mpicc as that would use GNU compilers!!!
export FC=$MKLROOT/../mpi/intel64/bin/mpiifort
export CC=$MKLROOT/../mpi/intel64/bin/mpiicc
export CXX=icpc

rm -r objects/* # one-time, if you build for Gfortran previously


# some systems don't have mpiifort for Intel
# use ifort as mpif90 get partially picked-up as GNU
OPTS="-DUSEMKL=yes -DLIB_DIR=$HOME/flibs-intel"

[[ $1 == "-d" ]] && OPTS="-DCMAKE_BUILD_TYPE=Debug $OPTS"

cmake $OPTS -B objects .

cmake --build objects -j

