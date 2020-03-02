cmake_minimum_required(VERSION 3.15)
# 3.12 for find_ROOT, 3.13 for a lot of user functions, 3.14 for check_fortran_source_runs & fetchcontent
# 3.15 for robust Python finding

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR "use cmake -B build or similar to avoid building in-source, which is messy")
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "default to Release")
endif()

project(gemini3d
  LANGUAGES C Fortran  # MUST include C language for Intel / MKL to work
  DESCRIPTION "3-D ionospheric model"
  HOMEPAGE_URL https://github.com/gemini3d/gemini
  VERSION 0.3.0)

enable_testing()

if(realbits EQUAL 32)
  message(STATUS "32-bit real precision")
  set(arith s)
else()
  message(STATUS "64-bit real precision")
  set(realbits 64)
  set(arith d)
endif()

option(autobuild "autobuild missing Lapack, Scalapack or Mumps" ON)
option(glow "use NCAR GLOW instead of Fang" ON)
option(matlab "enable Matlab tests, which are slow and duplicate fast Octave tests")
option(hdf5 "use HDF5 file I/O" on)
option(netcdf "use NetCDF file I/O" off)
option(metis "use METIS" OFF)
option(scotch "use Scotch" OFF)
option(openmp "Use OpenMP when building MUMPS" OFF)

if(netcdf AND hdf5)
  message(FATAL_ERROR "HDF5 and NetCDF together is a planned feature")
endif()

set(FETCHCONTENT_UPDATES_DISCONNECTED off)
# "on" allows live debugging from _deps directory
# "off" is default, always checks for updates

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules/)

include(${CMAKE_SOURCE_DIR}/cmake/compilers.cmake)
if(hdf5)
include(${CMAKE_SOURCE_DIR}/cmake/h5fortran.cmake)
elseif(netcdf)
include(${CMAKE_SOURCE_DIR}/cmake/nc4fortran.cmake)
endif()
include(${CMAKE_SOURCE_DIR}/cmake/utils.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/mumps.cmake)

# --- simple tests to help catch broken libs
add_subdirectory(tests)

# --- vendor libraries
add_subdirectory(vendor/msis00)

# --- GEMINI
add_subdirectory(ionization)
add_subdirectory(numerical)
add_subdirectory(io)
add_subdirectory(temporal)
add_subdirectory(neutral)

add_library(collision OBJECT collisions/collisions.f90)
target_link_libraries(collision PRIVATE const)

add_subdirectory(sources)

# --- other libs
add_library(multifluid OBJECT multifluid/multifluid.f90)
target_link_libraries(multifluid PRIVATE advec calculus collision const diffusion grid mesh ionization mpimod precipBCs sources timeutils)

# --- Main GEMINI executable
add_executable(gemini.bin gemini.f90)
target_link_libraries(gemini.bin PRIVATE advec calculus const collision diffusion fang grid io interp ionization mesh mpimod multifluid
  neutral pathlib PDEelliptic PDEparabolic potential precipBCs reader sources temporal timeutils
  ${LAPACK_LIBRARIES})  # necessary for self-compiled

#--------magnetic field calculation executable---------------
add_executable(magcalc.bin magcalc.f90)
target_link_libraries(magcalc.bin PRIVATE mpimod const grid mesh io pathlib reader timeutils)

# --- self-tests
include(${CMAKE_SOURCE_DIR}/cmake/test2d.cmake)
#if(NOT DEFINED ENV{CI} OR NOT $ENV{CI})
  # need the dollar sign on ENV for IF: save needless downloading
include(${CMAKE_SOURCE_DIR}/cmake/test3d.cmake)
#endif()

foreach(t collision fang fang_run io ionization logging multifluid neutral pathlib precipBCs reader sources temporal timeutils)
  target_include_directories(${t} INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
  set_target_properties(${t} PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
endforeach()