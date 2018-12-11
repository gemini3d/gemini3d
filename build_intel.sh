#!/bin/bash
# prereqs -- you'll need these in your current runtime shell as well to successfully run the program.
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
# ------ temporary environment follows
(

rm -r objects/* # one-time, if you build for Gfortran previously

cd objects

# some systems don't have mpiifort for Intel
# use ifort as mpif90 get partially picked-up as GNU
cmake \
-DUSEMKL=yes -DLIB_DIR=$HOME/flibs-intel .. \
#-DCMAKE_BUILD_TYPE=Debug

)

# Requires CMake 3.12
cmake --build objects -j

