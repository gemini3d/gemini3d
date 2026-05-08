# helps debug CMake version related compile or link issues by using the JSON Compilation Database output by Ninja.
# Note that GNU Make can only output compile commands, not link commands in the JSON Compilation Database at this time.
#
# this approach takes about 30 seconds per CMake version on a laptop

cmake_minimum_required(VERSION 3.25)

set(cmake_vers 3.25.3 3.26.5 3.27.9 3.28.6)

set(CMAKE_FIND_APPBUNDLE NEVER)
# must have this for macOS or it will open the CMake GUI and fail.


function(find_cmake cmake_req)

find_program(cmake_${cv}
NAMES cmake
PATHS ~/cmake-${cv}/
PATH_SUFFIXES bin CMake.app/Contents/bin
NO_DEFAULT_PATH
)

if(NOT cmake_${cv})
  return()
endif()

message(DEBUG "Found CMake at ${cmake_${cv}}")

execute_process(COMMAND ${cmake_${cv}} -E capabilities
OUTPUT_VARIABLE cmake_json
OUTPUT_STRIP_TRAILING_WHITESPACE
RESULT_VARIABLE _ret
)
if(NOT _ret EQUAL 0)
  message(WARNING "Failed to run ${cmake_${cv}} -E capabilities: ${cmake_json}")
  return()
endif()

message(TRACE "${cmake_json}")

string(JSON cmake_version ERROR_VARIABLE _err GET "${cmake_json}" "version" "string")
if(_err)
  message(WARNING "Failed to parse CMake version from JSON output of ${cmake_${cv}}: ${_err}")
  return()
endif()

if(NOT cmake_version STREQUAL "${cv}")
  message(WARNING "CMake version ${cv} found at ${cmake_${cv}}, but version output is ${cmake_version}. Skipping.")
  return()
endif()

message(STATUS "CMake version ${cv} found at ${cmake_${cv}}")

set(this_cmake ${cmake_${cv}} PARENT_SCOPE)

endfunction()


foreach(cv IN LISTS cmake_vers)

find_cmake(${cv})

set(target msis_ifc)

set(bindir ${CMAKE_CURRENT_LIST_DIR}/build-compdb)
# we use the same build dir to avoid all the paths being different in the compdb

execute_process(COMMAND ${this_cmake}
  -B ${bindir} -S ${CMAKE_CURRENT_LIST_DIR}/..
  -G Ninja
  -DCMAKE_BUILD_TYPE=Release
  --fresh
  RESULT_VARIABLE _ret
)
if(NOT _ret EQUAL 0)
  message(FATAL_ERROR "Failed to configure with CMake ${cv} at ${this_cmake}")
endif()

execute_process(COMMAND ${this_cmake} --build ${bindir} --verbose -- -t compdb-targets ${target}
  RESULT_VARIABLE _ret
  OUTPUT_VARIABLE compdb
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT _ret EQUAL 0)
  message(FATAL_ERROR "Failed to generate compdb with CMake ${cv} at ${this_cmake}")
  continue()
endif()

set(compdb_fn ${CMAKE_CURRENT_LIST_DIR}/compdb-${cv}.json)

message(STATUS "Writing compilation database for CMake ${cv} to ${compdb_fn}")

file(WRITE ${compdb_fn} "${compdb}")


endforeach()
