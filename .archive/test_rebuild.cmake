# script the testing of autobuild rebuilds
# these can develop problems with variables falling
# back to default on rebuild, with link-time failure

# let's do a basic test, as HDF5 will fail with this if link-time becomes broken

cmake_minimum_required(VERSION 3.20)

set(opts)
set(tgt msis_setup)

set(CTEST_USE_LAUNCHERS 1)
set(CTEST_TEST_TIMEOUT 10)
set(CTEST_OUTPUT_ON_FAILURE true)

set(_s ${CTEST_SCRIPT_DIRECTORY}/..)
cmake_path(ABSOLUTE_PATH _s NORMALIZE)
set(CTEST_SOURCE_DIRECTORY ${_s})
if(NOT DEFINED CTEST_BINARY_DIRECTORY)
  set(CTEST_BINARY_DIRECTORY ${CTEST_SOURCE_DIRECTORY}/build_test)
endif()

if(NOT CMAKE_BUILD_TYPE)
  # RelWithDebInfo -O2, Release -O3
  set(CMAKE_BUILD_TYPE Release)
endif()
list(APPEND opts -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})

include(${CTEST_SOURCE_DIRECTORY}/cmake/find_generator.cmake)
include(${CTEST_SOURCE_DIRECTORY}/cmake/cpu_count.cmake)

# -- run main test

ctest_start(Experimental)

ctest_configure(OPTIONS "${opts}"
RETURN_VALUE _ret
CAPTURE_CMAKE_ERROR _err)
if(NOT (_ret EQUAL 0 AND _err EQUAL 0))
  message(FATAL_ERROR "Build failed: return cod ${_ret} CMake Error ${_err}")
endif()

# ensure that rebuild works
ctest_build(TARGET ${tgt}
RETURN_VALUE _ret
CAPTURE_CMAKE_ERROR _err)
if(NOT (_ret EQUAL 0 AND _err EQUAL 0))
  message(FATAL_ERROR "Build failed: return cod ${_ret} CMake Error ${_err}")
endif()

ctest_build(TARGET ${tgt}
RETURN_VALUE _ret
CAPTURE_CMAKE_ERROR _err)
if(NOT (_ret EQUAL 0 AND _err EQUAL 0))
  message(FATAL_ERROR "RE-Build failed: return cod ${_ret} CMake Error ${_err}")
endif()
