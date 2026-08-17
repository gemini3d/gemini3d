#!/usr/bin/env bash
set -eo pipefail

gemrun=${1:-build/gemini3d.run}
datadir=$2

if [[ ! -d "$datadir" ]]; then
    echo "Data directory does not exist: $datadir"
    exit 1
fi

if [[ ! -f "$gemrun" ]]; then
    echo "gemini3d.run program does not exist: $gemrun"
    exit 1
fi

gemrun=$(realpath "$gemrun")
datadir=$(realpath "$datadir")
builddir=$(dirname "$gemrun")

case $OSTYPE in
  darwin*)  export DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH:$builddir/_deps/hdf5_zlib-build:$builddir/_deps/hdf5-build/bin" ;;
  linux*)   export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$builddir/_deps/hdf5_zlib-build:$builddir/_deps/hdf5-build/bin" ;;
  *)        echo "Unknown OS type: $OSTYPE"; exit 1 ;;
esac

# execute from build directory so the data files are found
(
cd "$builddir"

./gemini3d.run
./gemini3d.run "$datadir"
)
