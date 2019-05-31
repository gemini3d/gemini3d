#!/bin/bash

# Intel Fortran needs only Mumps.
# CMake >= 3.12 required

set -e  # abort on any error
set -u  # abort on undefined variable (a common bash goof)

# ======= user config ====================
# Arbitrary names/locations

# all libraries installed under $PREFIX/libraryname
PREFIX=$HOME/lib_gemini_intel18

# working directory, so you can rebuild later without complete recompilation
WD=/tmp

# true to erase/redo builds
wipe=false
# for each library, switch "true" to "false" if you don't want it.
BUILDMUMPS=true

# which compilers do you want?
export FC=$(which mpiifort)
export CC=$(which mpiicc)
export CXX=$(which icpc)

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
MUMPSGIT=https://github.com/scivision/mumps
MUMPSPREFIX=$PREFIX/mumps

cmake --version

mkdir -p $WD

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

cmake -DCMAKE_INSTALL_PREFIX=$MUMPSPREFIX -B $WD/mumps/build -S $WD/mumps

cmake --build $WD/mumps/build -j --target install -- $LOADLIMIT

fi
