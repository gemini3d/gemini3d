#!/bin/sh
# for HPC or elsewhere you need to compile everything except the compiler

set -u
set -e

R=$1
B=build

mpi_root=$R/lib_gcc/openmpi-3.1.5/bin

#=======================================
export FC=gfortran
export CC=gcc

cmake -B$B -DMPI_ROOT=$mpidir -Dhdf5=on

cmake --build $B -j
