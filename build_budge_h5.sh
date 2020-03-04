#!/bin/sh

set -e

builddir=build_budge_h5/

. $HOME/gcccompilers.sh

OPTS="-Dhdf5=on"

rm -rf $builddir

cmake -B $builddir $OPTS

cmake --build $builddir -j
