#!/bin/sh
# for HPC or elsewhere you need to compile everything except the compiler

set -u
set -e

R=$1
B=build

mpidir=$R/lib_gcc/openmpi-3.1.4/bin

#=======================================
mpidir=$(readlink -f $mpidir)

export FC=$mpidir/mpif90
export CC=$mpidir/mpicc
export MPIFC=$FC
export MPICC=$CC

meson $B -DMPI_ROOT=$mpidir -Drealbits=32 -Darith=s -Dsystem_blas=true

ninja -C $B
