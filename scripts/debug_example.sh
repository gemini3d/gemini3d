#!/bin/bash
# this is to debug simulations from gemini3d/gemini-examples repo
#
# Usage: ./debug_example.sh <config_nml_to_copy> <simulation_output_dir> [cmake_opts...]

[[ $# -lt 2 ]] && { echo "Usage: $0 <config_nml_to_copy> <simulation_output_dir> [cmake_opts...]"; exit 1; }

config_nml_to_copy=$1
simulation_output_dir=$2
if [[ $# -ge 3 ]]; then
  cmake_opts="${@:3}"
fi

[[ -f "$config_nml_to_copy" ]] || { echo "Config NML file '$config_nml_to_copy' does not exist"; exit 1; }

if [[ ! -d "$simulation_output_dir" ]]; then
  mkdir -p "$simulation_output_dir" || { echo "Failed to create simulation output directory '$simulation_output_dir'"; exit 1; }
fi

this_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

builddir="$simulation_output_dir/build"

copied_config="$simulation_output_dir/$(basename "$config_nml_to_copy")"
if [[ ! -f "$copied_config" ]] || ! diff -q -- "$config_nml_to_copy" "$copied_config" >/dev/null; then
  cp -v "$config_nml_to_copy" "$simulation_output_dir/" || { echo "Failed to copy '$config_nml_to_copy' to '$simulation_output_dir'"; exit 1; }
fi

cmake -S $this_script_dir/.. -B $builddir $cmake_opts

cmake --build $builddir
