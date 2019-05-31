#!/bin/bash

# this script is to help those who must compile everything except the compiler.
# For example, CentOS < 8 users.

# if using Gfortran, Gfortran >= 6 is REQUIRED
# CMake >= 3.12 required

# Newer Gfortran versions e.g. Gfortran 8, 9, etc. are recommended for better performance.
# Gfortran 6 is the oldest currently supported version, so stay ahead by using say Gfortran 8 or newer.


set -e  # abort on any error
set -u  # abort on undefined variable (a common bash goof)

# ======= user config ====================
# Arbitrary names/locations

# all libraries installed under $PREFIX/libraryname
PREFIX=$HOME/lib_gemini_gcc8

# working directory, so you can rebuild later without complete recompilation
WD=/tmp

# true to erase/redo builds
wipe=false
# for each library, switch "true" to "false" if you don't want it.
BUILDMPI=true
BUILDLAPACK=true
BUILDSCALAPACK=true
BUILDMUMPS=true

# which compilers do you want?
export FC=$(which gfortran)
export CC=$(which gcc)
export CXX=$(which g++)

# [optional] if you want to limit build load factor to this (for slow laptops)
# does NOT impact library performance, just how fast it builds.
# blank means maximum build speed
LOADLIMIT=
# LOADLIMIT="-l 4"

# ================================================
# normally don't adjust parameters below this line

echo "FC=$FC"
echo "CC=$CC"
echo "CXX=$CXX"

# Library parameters
LAPACKGIT=https://github.com/Reference-LAPACK/lapack
LAPACKPREFIX=$PREFIX/lapack

MPIVERSION=3.1.3  # MPI 4 doesn't currently work with ScalaPack
MPIFN=openmpi-$MPIVERSION.tar.bz2
MPIURL=https://download.open-mpi.org/release/open-mpi/v3.1/$MPIFN
MPISHA1=b3c60e2bdd5a8a8e758fd741f9a5bebb84da5e81
MPIPREFIX=$PREFIX/openmpi-$MPIVERSION

SCALAPACKGIT=https://github.com/scivision/scalapack
SCALAPACKPREFIX=$PREFIX/scalapack

MUMPSGIT=https://github.com/scivision/mumps
MUMPSPREFIX=$PREFIX/mumps

[[ $($FC -dumpversion) < 6 ]] && { echo "Gfortran >= 6 required"; exit 1; }

cmake --version

mkdir -p $WD

#==============================================================
# OpenMPI 3.1
# At this time, Scalapack isn't compatible with OpenMPI 4
# https://www.open-mpi.org/software/ompi/v3.1/

if $BUILDMPI
then

[[ -f $WD/$MPIFN ]] || curl -L $MPIURL -o $WD/$MPIFN

[[ $(sha1sum $WD/$MPIFN | cut -f1 -d' ') == $MPISHA1 ]] || { echo "checksum not match $WD/$MPIFN"; exit 1; }

tar -xf $WD/$MPIFN -C $WD

echo "installing OpenMPI $MPIVERSION to $MPIPREFIX"

(cd $WD/openmpi-$MPIVERSION; ./configure --prefix=$MPIPREFIX CC=$CC CXX=$CXX FC=$FC)

make -C $WD/openmpi-$MPIVERSION -j $LOADLIMIT

make -C $WD/openmpi-$MPIVERSION install

fi

#================================================================
# LAPACK

if $BUILDLAPACK
then

if [[ -d $WD/lapack ]]
then
git -C $WD/lapack pull
else
git clone --depth 1 $LAPACKGIT $WD/lapack
fi

mkdir $WD/lapack/build

if $wipe
then
rm -f $WD/lapack/build/CMakeCache.txt
rm -r $LAPACKPREFIX/*
fi

cmake -DCMAKE_INSTALL_LIBDIR=$LAPACKPREFIX -B $WD/lapack/build -S $WD/lapack
cmake --build $WD/lapack/build -j --target install -- $LOADLIMIT

fi

#===============
# scalapack

if $BUILDSCALAPACK
then

if [[ -d $WD/scalapack ]]
then
git -C $WD/scalapack pull
else
git clone --depth 1 $SCALAPACKGIT $WD/scalapack
fi

if $wipe
then
rm -f $WD/scalapack/build/CMakeCache.txt
rm -r $SCALAPACKPREFIX/*
fi

cmake -DCMAKE_INSTALL_PREFIX=$SCALAPACKPREFIX -DMPI_ROOT=$MPIPREFIX -DLAPACK_ROOT=$LAPACKPREFIX -B $WD/scalapack/build -S $WD/scalapack

cmake --build $WD/scalapack/build -j --target install -- $LOADLIMIT
fi

#=================
# MUMPS

if $BUILDMUMPS
then

if [[ -d $WD/mumps ]]
then
git -C $WD/mumps pull
else
git clone --depth 1 $MUMPSGIT $WD/mumps
fi

if $wipe
then
rm -f $WD/mumps/build/CMakeCache.txt
rm -r $MUMPSPREFIX/*
fi

cmake -DCMAKE_INSTALL_PREFIX=$MUMPSPREFIX -DSCALAPACK_ROOT=$SCALAPACKPREFIX -DMPI_ROOT=$MPIPREFIX -DLAPACK_ROOT=$LAPACKPREFIX -B $WD/mumps/build -S $WD/mumps

cmake --build $WD/mumps/build -j --target install -- $LOADLIMIT

fi
