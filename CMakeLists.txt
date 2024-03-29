cmake_minimum_required(VERSION 3.19...3.29)

# --- CMAKE_BUILD_TYPE default
# The simulations are 10x slower for default to Debug.
# Thus, for single config generators, set build type to Release
get_property(is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT is_multi_config AND NOT (CMAKE_BUILD_TYPE OR DEFINED ENV{CMAKE_BUILD_TYPE}))
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Release can be 10x faster simulation run time for gemini3d.run vs. Debug")
endif()

# --- auto toolchain for Cray
include(cmake/toolchain.cmake)

# --- main Gemini3D build

project(gemini3d
LANGUAGES C CXX Fortran
# Gemini3D is Fortran, but external libraries use C, and some find_package need C.
DESCRIPTION "3-D ionospheric model"
HOMEPAGE_URL https://github.com/gemini3d/gemini
VERSION 1.7.0
)

enable_testing()  # keep this so BUILD_TESTING=off doesn't remove all tests
include(CTest)

set(CMAKE_CXX_STANDARD 20)

include(cmake/package/git_rev.cmake)
include(cmake/GitSubmodule.cmake)
include(options.cmake)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Fortran filesystem library "ffilesystem"
git_submodule(${CMAKE_CURRENT_SOURCE_DIR}/ffilesystem)
report_submodule(${CMAKE_CURRENT_SOURCE_DIR}/ffilesystem)
add_subdirectory(ffilesystem)

# find python before excluding Anaconda
include(cmake/python.cmake)

# --- MPI
# NOTE: our find_package(MPI) needs to be before find_package(MUMPS), which also calls find_package(MPI)
# Conda e.g. scikit-learn can cause problems with finding MPI, so exclude Conda from search
# --- avoid Anaconda libraries
if(DEFINED ENV{CONDA_PREFIX})
  list(APPEND CMAKE_IGNORE_PREFIX_PATH $ENV{CONDA_PREFIX})
  list(APPEND CMAKE_IGNORE_PATH $ENV{CONDA_PREFIX}/bin)
  # need CMAKE_IGNORE_PATH for CMake < 3.23
  # and to ensure system env var PATH doesn't interfere
  # despite CMAKE_IGNORE_PREFIX_PATH
endif()

include(cmake/mpi.cmake)

# --- end MPI

git_submodule(${CMAKE_CURRENT_SOURCE_DIR}/h5fortran)
report_submodule(${CMAKE_CURRENT_SOURCE_DIR}/h5fortran)
add_subdirectory(h5fortran)
message(STATUS "${PROJECT_NAME}: HDF5 ${HDF5_VERSION}")
# this also implicitly finds HDF5.
# h5fortran is a high-level, object-oriented HDF5 interface.

# HDF5 bug #3663 for HDF5 1.14.2, 1.14.3, ...?
# https://github.com/HDFGroup/hdf5/issues/3663
if(WIN32 AND CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
if(HDF5_VERSION MATCHES "1.14.[2-3]")
  message(VERBOSE "HDF5: applying workaround for HDF5 bug #3663 with Intel oneAPI on Windows")
  set(oneapi_hdf5_bug3663 shlwapi)
endif()
endif()

include(cmake/compilers.cmake)

# --- linear algebra libraries
if(NOT DEFINED ENV{MKLROOT})
  find_package(LAPACK REQUIRED)
endif()
find_package(SCALAPACK REQUIRED)

find_package(MUMPS CONFIG)
if(NOT MUMPS_FOUND)
  message(FATAL_ERROR ${need_gemext})
endif()

# --- climate models
if(glow)
  git_submodule(${CMAKE_CURRENT_SOURCE_DIR}/glow)
  report_submodule(${CMAKE_CURRENT_SOURCE_DIR}/glow)
  add_subdirectory(glow)
endif()

if(hwm14)
  git_submodule(${CMAKE_CURRENT_SOURCE_DIR}/hwm14)
  report_submodule(${CMAKE_CURRENT_SOURCE_DIR}/hwm14)
  add_subdirectory(hwm14)
endif()

# --- MSISE00 / MSIS 2.x
set(msis2 on)
set(msis_BUILD_UTILS on)
set(msis_BUILD_TESTING on) # so we know MSIS is working as it's fundamental to Gemini3D
git_submodule(${CMAKE_CURRENT_SOURCE_DIR}/msis)
report_submodule(${CMAKE_CURRENT_SOURCE_DIR}/msis)
add_subdirectory(msis)

# --- Gemini3D library
add_subdirectory(src)

# --- Gemini3D self test
include(cmake/test/launcher.cmake) # for MPI tests

add_subdirectory(test)
# fundamental tests of MPI and numeric libraries essential for Gemini3D

# self-test simulations -- after all targets for if(TARGET ...)
include(cmake/test/config.cmake)
include(cmake/test/sim.cmake)

# summary print
include(cmake/summary.cmake)

# packaging
include(cmake/package/pkgconf.cmake)
include(cmake/package/install.cmake)
