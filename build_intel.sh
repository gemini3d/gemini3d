#!/bin/bash
# prereqs -- you'll need these in your current runtime shell as well to successfully run the program.
# *** for subsequent builds, you can just type "make" in the objects/ directory ***
# (I keep a second Terminal tab for this purpose)
# I.e. copy and paste them or put in ~/.bashrc
#----- gemini -------------
(
. $MKLROOT/../bin/compilervars.sh intel64
. $MKLROOT/bin/mklvars.sh intel64 ilp64

rm -r objects/* # one-time, if you build for Gfortran previously
cd objects

# some systems don't have mpiifort for Intel
# use ifort as mpif90 get partially picked-up as GNU
FC=$MKLROOT/../mpi/intel64/bin/mpiifort CC=$MKLROOT/../mpi/intel64/bin/mpiicc cmake -DUSEMKL=yes -DLIB_DIR=~/flibs-intel ..

)

# Requires CMake 3.12
cmake --build objects -j

