#!/bin/sh

set -e

builddir=build_budge/

. $HOME/gcccompilers.sh

OPTS="-Dhdf5=off"

cmake -B $builddir $OPTS

cmake --build $builddir -j
