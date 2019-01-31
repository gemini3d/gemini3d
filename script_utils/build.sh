#!/bin/bash

# helper script, not to run directly

set -u
set -e

WD=$(dirname ${BASH_SOURCE[0]})

#. $WD/meson.sh

#[[ $MESONSTAT == 0 ]] || 

. $WD/cmake.sh


