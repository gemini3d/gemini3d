#!/bin/bash

# this is helper script, not meant to be run directly

set -u
set -e

cmake --version

[[ ${1:-} == "-d" ]] && OPTS="-DCMAKE_BUILD_TYPE=Debug $OPTS"
[[ ${1:-} == "-t" ]] && OPTS="-DTRACE:BOOL=on $OPTS"

cmake $OPTS -B objects -S .

cmake --build objects -j
