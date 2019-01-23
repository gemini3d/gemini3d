#!/bin/bash

# this script is to help those who must compile everything except the compiler.
# For example, CentOS < 8 users.

# if using Gfortran, Gfortran >= 6 is REQUIRED (newer versions e.g. Gfortran 7, 8, 9, etc. are recommended for better performance)
#
# Gemini currently has a bug with Ifort, but once that's fixed, Ifort should work. For now use Gfortran.
#
# Flang, PGI and/or NAG support are anticipated soon from vendors, possibly in 2019. 
# Ask if desired.

#
# for each library, switch "true" to "false" if you don't want it.

set -e  # abort on any error
set -u  # abort on undefined variable (a common bash goof)

export FC=gfortran
export CC=gcc
export CXX=g++

WD=$HOME/libs_gemini  # so you can rebuild later without complete recompilation
PREFIX=$HOME/.local  # all libraries installed under this directory

[[ $(gfortran -dumpversion) < 6 ]] && { echo "Gfortran >= 6 required"; exit 1; }

mkdir -p $WD

#==============================================================
# OpenMPI 3.1
# At this time, Scalapack isn't compatible with OpenMPI 4
# https://www.open-mpi.org/software/ompi/v3.1/
MPIVERSION=3.1.3
MPIFN=openmpi-$MPIVERSION.tar.bz2
MPIURL=https://download.open-mpi.org/release/open-mpi/v3.1/$MPIFN
MPISHA1=b3c60e2bdd5a8a8e758fd741f9a5bebb84da5e81
MPIPREFIX=$PREFIX/openmpi-$MPIVERSION-gcc$(gfortran -dumpversion)

if true; then

cd $WD

[[ -f $MPIFN ]] || curl -L $MPIURL -o $MPIFN

[[ $(sha1sum $MPIFN | cut -f1 -d' ') == $MPISHA1 ]] || { echo "checksum not match $MPIFN"; exit 1; }

tar -xf $MPIFN
cd $WD/openmpi-$MPIVERSION

echo "installing OpenMPI $MPIVERSION to $MPIPREFIX"

./configure --prefix=$MPIPREFIX CC=gcc CXX=g++ FC=gfortran

make -j -l2

make install

fi

#================================================================
# LAPACK
LAPACKGIT=https://github.com/Reference-LAPACK/lapack
LAPACKPREFIX=$PREFIX/lapack-gcc$(gfortran -dumpversion)

if true; then

cd $WD

git clone --depth 1 $LAPACKGIT

cd lapack
mkdir build
cd build
cmake -DCMAKE_INSTALL_LIBDIR=$LAPACKPREFIX ..
cmake --build -j . --target install -- -l 2

fi

#===============
# scalapack

MUMPSGIT=https://github.com/scivision/fortran-libs

if true; then

cd $WD

git clone --depth 1 $MUMPSGIT mumps

cd mumps

./build_self.sh

fi
