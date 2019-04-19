#!/bin/bash

# this is helper script, not meant to be run directly

set -u
set -e

MPIPREFIX=$PREFIX/openmpi-3.1.3
LAPACKPREFIX=$PREFIX/lapack
SCALAPACKPREFIX=$PREFIX/scalapack
MUMPSPREFIX=$PREFIX/mumps

export FC=$MPIPREFIX/bin/mpifort
export CC=$MPIPREFIX/bin/mpicc

# ============================================================
# this temporarily disables Intel compiler (if installed) from messing up your gfortran environment.
MKLROOT=
LD_LIBRARY_PATH=
