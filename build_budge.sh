#!/bin/sh

set -e

. $HOME/gcccompilers.sh

OPTS="-DMPI_ROOT=$HOME/projs/lib_gcc/openmpi-3.1.5 -Dhdf5=off"

rm -rf build_budge/

cmake -B build_budge $OPTS

cmake --build build -j