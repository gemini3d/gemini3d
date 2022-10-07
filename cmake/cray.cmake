# toolchain for Intel oneAPI and/or GCC compilers on Cray system
# canNOT use from Project CMakeLists.txt.
# To propagate to ExternalProject, cannot use any "-D" variables or set(ENV{}) in this file.
#
# NOTE: your Cray system may have different versions/paths, treat this like a template.
#
# Copy this file to a convenient location like ~ directory, and use with any project on Cray like:
#  cmake -DCMAKE_TOOLCHAIN_FILE=~/cray.cmake -B build

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
if(NOT ret EQUAL 0)
  message(WARNING "failed to load ${gcc_mod}:    ${out}  ${ret}")
endif()

find_program(cc NAMES gcc REQUIRED)

execute_process(COMMAND ${cc} -dumpversion
OUTPUT_VARIABLE gcc_vers
ERROR_VARIABLE err
RESULT_VARIABLE ret
)
if(NOT ret EQUAL 0)
  message(WARNING "ERROR: failed to get ${gcc_mod} version: ${ret}  ${err}")
  return()
endif()
if(gcc_vers VERSION_LESS 9.1)
  message(WARNING "GCC toolchain >= 9.1 is required for oneAPI")
  return()
endif()

execute_process(COMMAND ${cc} -v
OUTPUT_VARIABLE gcc_verb
ERROR_VARIABLE gcc_verb
RESULT_VARIABLE ret
)
if(NOT ret EQUAL 0)
  message(WARNING "ERROR: failed to get ${gcc_mod} build details: ${ret}   ${gcc_verb}")
  return()
endif()

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
  message(WARNING "Unknown toolchain program environment on ${host}.
  Check that module names are set for this system in ${CMAKE_CURRENT_LIST_FILE}:
  ${pegnu} ${peintel} ${pecray}")
endif()
