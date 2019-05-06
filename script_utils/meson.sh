#!/bin/bash

# this is helper script, not meant to be run directly

meson --version

[[ ${1:-} == "-d" ]] && OPTS="--buildtype=debug $OPTS" || OPTS="--buildtype=release $OPTS"
[[ ${1:-} == "-t" ]] && OPTS="-DTRACE=on $OPTS"


meson $OPTS setup objects

MESONSTAT=$?

ninja -C objects
