#!/bin/bash

# this script is to help those who must compile everything except the compiler.
# For example, CentOS < 8 users.

# if using Gfortran, Gfortran >= 6 is REQUIRED
# Newer Gfortran versions e.g. Gfortran 8, 9, etc. are recommended for better performance.
# Gfortran 6 is the oldest currently supported version, so stay ahead by using say Gfortran 8 or newer.
#
# Gemini currently has a bug with Ifort, but once that's fixed, Ifort should work. For now use Gfortran.
#
# Flang, PGI and/or NAG support are anticipated soon from vendors, possibly in 2019.
# Ask if desired.

set -e  # abort on any error
set -u  # abort on undefined variable (a common bash goof)

# ======= user config ====================
# Arbitrary names/locations

# all libraries installed under $PREFIX/libraryname
PREFIX=$HOME/.local
# whatever name you want to name at end of each library directory, arbitrary
SUFFIX=gcc8
# working directory, so you can rebuild later without complete recompilation
WD=$HOME/libs_gemini-$SUFFIX

# for each library, switch "true" to "false" if you don't want it.
BUILDMPI=false
BUILDLAPACK=false
BUILDSCALAPACK=true
BUILDMUMPS=true

# which compilers do you want?
export FC=$(which gfortran)
export CC=$(which gcc)
export CXX=$(which g++)

# ================================================
# normally don't adjust parameters below this line

echo "FC=$FC"
echo "CC=$CC"
echo "CXX=$CXX"

# Library parameters
LAPACKGIT=https://github.com/Reference-LAPACK/lapack
LAPACKPREFIX=$PREFIX/lapack-$SUFFIX

MPIVERSION=3.1.3  # MPI 4 doesn't currently work with ScalaPack
MPIFN=openmpi-$MPIVERSION.tar.bz2
MPIURL=https://download.open-mpi.org/release/open-mpi/v3.1/$MPIFN
MPISHA1=b3c60e2bdd5a8a8e758fd741f9a5bebb84da5e81
MPIPREFIX=$PREFIX/openmpi-$MPIVERSION-$SUFFIX

SCALAPACKGIT=https://github.com/scivision/scalapack
SCALAPACKPREFIX=$PREFIX/scalapack-$SUFFIX

MUMPSGIT=https://github.com/scivision/mumps
MUMPSPREFIX=$PREFIX/mumps-$SUFFIX

[[ $($FC -dumpversion) < 6 ]] && { echo "Gfortran >= 6 required"; exit 1; }

cmake --version

mkdir -p $WD

#==============================================================
# OpenMPI 3.1
# At this time, Scalapack isn't compatible with OpenMPI 4
# https://www.open-mpi.org/software/ompi/v3.1/

if $BUILDMPI
then

cd $WD

[[ -f $MPIFN ]] || curl -L $MPIURL -o $MPIFN

[[ $(sha1sum $MPIFN | cut -f1 -d' ') == $MPISHA1 ]] || { echo "checksum not match $MPIFN"; exit 1; }

tar -xf $MPIFN
cd $WD/openmpi-$MPIVERSION

echo "installing OpenMPI $MPIVERSION to $MPIPREFIX"

./configure --prefix=$MPIPREFIX CC=$CC CXX=$CXX FC=$FC

make -j -l 4

make install

fi

#================================================================
# LAPACK

if $BUILDLAPACK
then

cd $WD

[[ -d lapack ]] && { (cd lapack; git pull) } || git clone --depth 1 $LAPACKGIT

cd lapack
mkdir -p build

cmake -DCMAKE_INSTALL_LIBDIR=$LAPACKPREFIX -B build -S .
cmake --build build -j --target install -- -l 4

fi

#===============
# scalapack

if $BUILDSCALAPACK
then

cd $WD

[[ -d scalapack ]] && { (cd scalapack; git pull) } || git clone --depth 1 $SCALAPACKGIT

cmake -DCMAKE_INSTALL_PREFIX=$SCALAPACKPREFIX -DMPI_ROOT=$MPIPREFIX -DLAPACK_ROOT=$LAPACKPREFIX -B scalapack/build -S scalapack

cmake --build scalapack/build -j --target install -- -l 4

fi

#=================
# MUMPS

if $BUILDMUMPS; then

cd $WD

[[ -d mumps ]] && { (cd mumps; git pull) } || git clone --depth 1 $MUMPSGIT mumps

cmake -DCMAKE_INSTALL_PREFIX=$MUMPSPREFIX -DSCALAPACK_ROOT=$SCALAPACKPREFIX -DMPI_ROOT=$MPIPREFIX -DLAPACK_ROOT=$LAPACKPREFIX -B mumps/build -S mumps

cmake --build mumps/build -j --target install -- -l 4

fi
