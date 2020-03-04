#!/bin/sh

set -e

builddir=build_budge/

. $HOME/gcccompilers.sh

OPTS="-Dhdf5=off"

rm -rf $builddir

cmake -B $builddir $OPTS

cmake --build $builddir -j
