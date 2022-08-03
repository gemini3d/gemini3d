# loads modules on Cray system
# canNOT use from Project CMakeLists.txt.
# To propagate to ExternalProject, cannot use any "-D" variables or set(ENV{}) in this file.
#
# NOTE: your Cray system may have different versions/paths, treat this like a template.
#
# NOTE: to specify install directory, run like:
#   cmake -DCMAKE_INSTALL_PREFIX=<install_dir> -P cray.cmake

cmake_minimum_required(VERSION 3.20.2)

# --- module names (may be different on your system)

set(gcc_mod gcc)
set(pecray PrgEnv-cray)
set(pegnu PrgEnv-gnu)
set(peintel PrgEnv-intel)

# --- main script

find_package(EnvModules REQUIRED)

function(gcc_toolchain)

set(CXXFLAGS $ENV{CXXFLAGS})
if(CXXFLAGS MATCHES "--gcc-toolchain")
  return()
endif()

env_module(load ${gcc_mod} OUTPUT_VARIABLE out RESULT_VARIABLE ret)
message(STATUS "load ${gcc_mod}:    ${out}  ${ret}")

find_program(cc NAMES gcc REQUIRED)

execute_process(COMMAND ${cc} -dumpversion
OUTPUT_VARIABLE gcc_vers
COMMAND_ERROR_IS_FATAL ANY
)
if(gcc_vers VERSION_LESS 9.1)
  message(FATAL_ERROR "GCC toolchain >= 9.1 is required for oneAPI")
endif()

execute_process(COMMAND ${cc} -v
OUTPUT_VARIABLE gcc_verb
ERROR_VARIABLE gcc_verb
COMMAND_ERROR_IS_FATAL ANY
)

set(pat "--prefix=([/a-zA-Z0-9_\\-\\.]+)")
string(REGEX MATCH "${pat}" gcc_prefix "${gcc_verb}")

if(CMAKE_MATCH_1)
  string(APPEND CXXFLAGS " --gcc-toolchain=${CMAKE_MATCH_1}")
  set(ENV{CXXFLAGS} ${CXXFLAGS})
else()
  message(WARNING "GCC toolchain not found")
endif()

endfunction(gcc_toolchain)

# the module commands only affect the current process, not the parent shell
env_module_list(mods)

cmake_host_system_information(RESULT host QUERY HOSTNAME)

if(mods MATCHES "${peintel}")
  # check GCC toolchain for Intel compiler
  gcc_toolchain()
  message(STATUS "Using ${peintel} program environment on ${host}")
elseif(mods MATCHES "${pegnu}")
  message(STATUS "Using ${pegnu} program environment on ${host}")
elseif(mods MATCHES "${pecray}")
  message(WARNING "${pecray} PE may not work with this project on ${host}. Try command like this first:
  module swap ${pecray} ${pegnu}
  OR
  module swap ${pecray} ${peintel}")
else()
  message(WARNING "Unknown toolchain program environment on ${host}")
endif()
