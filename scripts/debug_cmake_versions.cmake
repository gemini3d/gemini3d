# helps debug CMake version related compile or link issues by using the JSON Compilation Database output by Ninja.
# Note that GNU Make can only output compile commands, not link commands in the JSON Compilation Database at this time.
#
# this approach takes about 30 seconds per CMake version on a laptop

cmake_minimum_required(VERSION 3.26)

set(cmake_vers 3.26.6 3.27.9 3.28.6)

set(cmake_gen Ninja)
if(WIN32)
  list(APPEND cmake_gen "MinGW Makefiles")
else()
  list(APPEND cmake_gen "Unix Makefiles")
endif()

set(CMAKE_FIND_APPBUNDLE NEVER)
# must have this for macOS or it will open the CMake GUI and fail.

set(compdb_target msis_ifc)
# this is a target that the build will fail on for CMake < 3.28

set(mod_remove _deps/msis-build/include/msis_constants.mod)


function(find_cmake cmake_req cmake_var)

find_program(cmake_${cv}
NAMES cmake
PATHS ~/cmake-${cv}/
PATH_SUFFIXES bin CMake.app/Contents/bin
REQUIRED
NO_DEFAULT_PATH
)

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

set(${cmake_var} ${cmake_${cv}} PARENT_SCOPE)

endfunction()

# Main Loop

set(cmake_failed_versions)
foreach(cv IN LISTS cmake_vers)
  foreach(gen IN LISTS cmake_gen)

find_cmake(${cv} this_cmake)

cmake_host_system_information(RESULT Ncpu QUERY NUMBER_OF_PHYSICAL_CORES)

set(srcdir ${CMAKE_CURRENT_LIST_DIR}/..)
set(bindir ${srcdir}/build-compdb-${gen})
# we use the same build dir to avoid all the paths being different in the compdb

execute_process(COMMAND ${this_cmake}
  -B ${bindir} -S ${srcdir}
  -G ${gen}
  -Dgemini3d_msis2:BOOL=on
  -DCMAKE_BUILD_TYPE=Release
  --fresh
  COMMAND_ECHO STDOUT
  RESULT_VARIABLE _ret
)
# gemini3d_msis2=on needed to trigger failure with older CMake, due to build graph generation incorrectness
if(NOT _ret EQUAL 0)
  message(FATAL_ERROR "Failed to configure with CMake ${cv} at ${this_cmake} with ${gen}
  ${_ret}")
endif()

# Write the Compilation Database for the specified target
if(gen STREQUAL "Ninja")
  execute_process(COMMAND ${this_cmake} --build ${bindir} --verbose -- -t compdb-targets ${compdb_target}
    RESULT_VARIABLE _ret
    OUTPUT_VARIABLE compdb
    COMMAND_ECHO STDOUT
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(NOT _ret EQUAL 0)
    message(FATAL_ERROR "Failed to generate compdb with CMake ${cv} at ${this_cmake}")
  endif()

  set(compdb_fn ${CMAKE_CURRENT_LIST_DIR}/compdb-${cv}.json)
  message(STATUS "Writing compilation database for CMake ${cv} to ${compdb_fn}")
  file(WRITE ${compdb_fn} "${compdb}")
endif()

# Attempt to build - expected to fail for CMake < 3.28 on Gemini3D - doesn't always fail due to random build order,
# so try a few times for failure

set(Nbuild 5) # arbitrary number of attempts to build

set(_par)
if(gen MATCHES "Makefiles")
  set(_par --parallel ${Ncpu})
endif()

foreach(i RANGE ${Nbuild})
  file(REMOVE ${bindir}/${mod_remove})

  execute_process(COMMAND ${this_cmake} --build ${bindir} ${_par} --target ${compdb_target} --clean-first
    RESULT_VARIABLE _ret
    COMMAND_ECHO STDOUT
  )
  if(_ret EQUAL 0)
    message(STATUS "SUCCESS: Built ${compdb_target} with CMake ${cv} at ${this_cmake} on attempt ${i} with ${gen}")
  else()
    message(WARNING "FAILURE: failed to build ${compdb_target} with CMake ${cv} at ${this_cmake} on attempt ${i} with ${gen}")
    list(APPEND cmake_failed_versions ${cv}-${gen})
    break()
  endif()
endforeach()

endforeach()
endforeach()


list(LENGTH cmake_failed_versions fail_count)
if(fail_count EQUAL 0)
  message(STATUS "SUCCESS: built ${compdb_target} with all CMake versions tested: ${cmake_vers}")
else()
  message(WARNING "FAILURE: CMake build failed with CMake versions: ${cmake_failed_versions} out of versions tested: ${cmake_vers}")
endif()
