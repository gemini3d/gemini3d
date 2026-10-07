#!/bin/bash
# this is to debug simulations e.g. from gemini3d/gemini-examples repo
#
# Usage: ./debug_example.sh <config_nml_to_copy> <simulation_output_dir> [cmake_opts...]
#
# Example from gemini3d/ directory:
#
#   scripts/debug_example.sh ../gemini-examples/init/CI-staging/varHWM/eq/config.nml $HOME/gemci/varHWM


[[ $# -lt 2 ]] && { echo "Usage: $0 <config_nml_to_copy> <simulation_output_dir> [cmake_opts...]"; exit 1; }

config_nml_to_copy=$1
simulation_output_dir=$2
if [[ $# -ge 3 ]]; then
  cmake_opts="${@:3}"
fi

# if CMAKE_BUILD_TYPE not set in variable cmake_opts, default to Release
if [[ ! "$cmake_opts" =~ -DCMAKE_BUILD_TYPE= ]]; then
  cmake_opts="-DCMAKE_BUILD_TYPE=Release $cmake_opts"
fi

[[ -f "$config_nml_to_copy" ]] || { echo "Config NML file '$config_nml_to_copy' does not exist"; exit 1; }

if [[ ! -d "$simulation_output_dir" ]]; then
  mkdir -p "$simulation_output_dir" || { echo "Failed to create simulation output directory '$simulation_output_dir'"; exit 1; }
fi

this_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

builddir="$simulation_output_dir/build"

copied_config="$simulation_output_dir/$(basename "$config_nml_to_copy")"
if [[ ! -f "$copied_config" ]] || ! diff -q -- "$config_nml_to_copy" "$copied_config" >/dev/null; then
  \cp -v "$config_nml_to_copy" "$simulation_output_dir/" || { echo "Failed to copy '$config_nml_to_copy' to '$simulation_output_dir'"; exit 1; }
fi

if [[ ! -f "$builddir/gemini.bin" ]]; then

cmake=$(command -v cmake) || { echo "CMake not found or working" >&2; exit 1; }

$cmake -S "$this_script_dir/.." -B "$builddir" $cmake_opts ||
  { echo "CMake configuration failed" >&2; exit 1; }

$cmake --build "$builddir" ||
  { echo "CMake build failed" >&2; exit 1; }

fi

case $OSTYPE in
  darwin*)
    Ncpu=$(sysctl -n hw.physicalcpu)
    ;;
  linux*)
    Ncpu=$(lscpu -p=CORE,SOCKET |
    awk -F, '!/^#/ { cores[$1 "," $2]=1 }
            END { for (core in cores) n++; print n }')
    ;;
  *)
    echo "Cannot determine the number of CPUs on: $OSTYPE; using 1" >&2
    Ncpu=1
    ;;
esac

[[ $Ncpu =~ ^[1-9][0-9]*$ ]] || {
  echo "Failed to determine physical CPU count" >&2
  exit 1
}

if [[ ! -f $simulation_output_dir/inputs/initial_conditions.h5 ]]; then
  py3=$(command -v python3) ||
    { echo "Python not found" >&2; exit 1; }

  "$py3" -m gemini3d --help ||
    { echo "PyGemini not available - get it from https://github.com/gemini3d/pygemini" >&2; exit 1; }

  echo "Setting up simulation inputs in '$simulation_output_dir'..."

  GEMINI_ROOT="$builddir" "$py3" -m gemini3d.model "$copied_config" "$simulation_output_dir" ||
    { echo "Failed to set up simulation inputs in '$simulation_output_dir'" >&2; exit 1; }
fi

echo "running simulation with $Ncpu MPI workers..."

mpiexec -np $Ncpu -wdir "$builddir" "$builddir/gemini.bin" "$simulation_output_dir"
