#!/bin/bash

set -u
set -e

MPIPREFIX=$MKLROOT/../mpi/intel64
LAPACKPREFIX=
SCALAPACKPREFIX=
MUMPSPREFIX=$PREFIX/mumps


export FC=$MPIPREFIX/bin/$FC
export CC=$MPIPREFIX/bin/$CC
export CXX=icpc
