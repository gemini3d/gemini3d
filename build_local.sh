#!/bin/sh
# for HPC or elsewhere you need to compile everything except the compiler

set -u
set -e

R=$1
B=build

mpidir=$R/lib_gcc/openmpi-3.1.4/bin

#=======================================
export FC=$mpidir/mpif90
export CC=$mpidir/mpicc
export MPIFC=$FC
export MPICC=$CC

meson $B -DMPI_ROOT=$mpidir

ninja -C $B
