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
VERSION 1.6.1
)

enable_testing()  # keep this so BUILD_TESTING=off doesn't remove all tests
include(CTest)

if(NOT DEFINED ${PROJECT_NAME}_BUILD_TESTING)
  set(${PROJECT_NAME}_BUILD_TESTING ${BUILD_TESTING})
endif()

include(cmake/package/git_rev.cmake)
include(options.cmake)

# Fortran filesystem library "ffilesystem"
find_package(ffilesystem CONFIG)
# this is the first package to find, so make a special error message for users
if(NOT ffilesystem_FOUND)
  message(FATAL_ERROR ${need_gemext})
endif()

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

find_package(h5fortran CONFIG REQUIRED)
# this also implicitly finds HDF5.
# h5fortran is a high-level, object-oriented HDF5 interface.

include(cmake/compilers.cmake)

# --- linear algebra libraries
if(NOT DEFINED ENV{MKLROOT})
  find_package(LAPACK REQUIRED)
endif()
find_package(SCALAPACK REQUIRED)
find_package(MUMPS REQUIRED CONFIG)

# --- climate models
if(glow)
  find_package(glow CONFIG REQUIRED)
endif()

if(hwm14)
  find_package(hwm14 CONFIG REQUIRED)
endif()

# --- MSISE00 or MSIS 2.x
find_package(msis CONFIG REQUIRED)
# WORKAROUND: MSIS shared Rpath fail
# Rpath is defined in MSIS, but ld fails to use it
get_target_property(msis_loc msis::msis_ifc LOCATION)
get_filename_component(msis_dir "${msis_loc}" DIRECTORY)

# END WORKAROUND

# unit test DLL paths
include(cmake/DllTestPath.cmake)


add_subdirectory(src)
# Gemini3D source code

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
