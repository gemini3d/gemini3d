cmake_minimum_required(VERSION 3.15)
# 3.14 for check_fortran_source_runs & fetchcontent
# 3.15 for robust Python finding

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR "use cmake -B build or similar to avoid building in-source, which is messy")
endif()

set(CMAKE_CONFIGURATION_TYPES "Release;RelWithDebInfo;Debug" CACHE STRING "Build type selections" FORCE)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "default build type")
endif()

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.18)
  cmake_policy(SET CMP0107 NEW)
  cmake_policy(SET CMP0108 NEW)
endif()
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.17)
  cmake_policy(SET CMP0099 NEW)
endif()

project(gemini3d
  LANGUAGES C Fortran  # MUST include C language for Intel / MKL to work
  DESCRIPTION "3-D ionospheric model"
  HOMEPAGE_URL https://github.com/gemini3d/gemini
  VERSION 0.6.1)


include(CTest) # for CDash
if(NOT DEFINED ${PROJECT_NAME}_BUILD_TESTING)
  set(${PROJECT_NAME}_BUILD_TESTING ${BUILD_TESTING})
endif()

if(realbits EQUAL 32)
  message(VERBOSE " 32-bit real precision")
  set(arith s)
else()
  message(VERBOSE " 64-bit real precision")
  set(realbits 64)
  set(arith d)
endif()

option(mpi "Use MPI parallelization" on)

option(autobuild "autobuild missing Lapack, Scalapack or Mumps" on)
option(glow "use NCAR GLOW airglow / aurora model" on)
option(hdf5 "use HDF5 file I/O" on)
option(netcdf "use NetCDF file I/O" off)
# MUMPS build options (only used if auto-building MUMPS)
option(metis "MUMPS: use METIS" off)
option(scotch "MUMPS: use Scotch" off)
option(openmp "MUMPS: use OpenMP" off)

# on: debug, off: normal
set(FETCHCONTENT_UPDATES_DISCONNECTED off)

# this helps linters e.g. Visual Studio Intellicode work properly
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/python.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/h5fortran.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/nc4fortran.cmake)

# do these after h5fortran to avoid h5fortran HDF5OK threads test failure
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/compilers.cmake)
# first this to get libs
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/mumps.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/compiler_flags.cmake)
# this comes last to avoid intermittent failures

# self-test simulations
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/test_setup.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/test_sim.cmake)

add_subdirectory(src)

# summary print
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/summary.cmake)
