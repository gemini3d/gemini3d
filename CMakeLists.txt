cmake_minimum_required(VERSION 3.14)
# 3.12 for find_ROOT, 3.13 for a lot of user functions, 3.14 for check_fortran_source_runs & fetchcontent
# 3.15 for robust Python finding

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR "use cmake -B build or similar to avoid building in-source, which is messy")
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "default to Release")
endif()

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.15)
  cmake_policy(SET CMP0094 NEW)
endif()

project(gemini3d
  LANGUAGES C Fortran  # MUST include C language for Intel / MKL to work
  DESCRIPTION "3-D ionospheric model"
  HOMEPAGE_URL https://github.com/gemini3d/gemini
  VERSION 0.3.3)

enable_testing()

if(realbits EQUAL 32)
  message(VERBOSE " 32-bit real precision")
  set(arith s)
else()
  message(VERBOSE " 64-bit real precision")
  set(realbits 64)
  set(arith d)
endif()

option(autobuild "autobuild missing Lapack, Scalapack or Mumps" on)
option(glow "use NCAR GLOW airglow / aurora model" on)
option(matlab "enable Matlab self-tests")
option(hdf5 "use HDF5 file I/O" on)
option(netcdf "use NetCDF file I/O" off)
# MUMPS build options (only used if auto-building MUMPS)
option(metis "MUMPS: use METIS" off)
option(scotch "MUMPS: use Scotch" off)
option(openmp "MUMPS: use OpenMP" off)

if(netcdf AND hdf5)
  message(FATAL_ERROR "HDF5 and NetCDF together is a planned feature")
endif()

# on: debug, off: normal
set(FETCHCONTENT_UPDATES_DISCONNECTED off)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/compilers.cmake)
if(hdf5)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/h5fortran.cmake)
elseif(netcdf)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/nc4fortran.cmake)
endif()
# do this last, after h5fortran to avoid Threads::Threads test failure
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/mumps.cmake)

# self-test simulations
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/test_setup.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/test_sim.cmake)

add_subdirectory(src)