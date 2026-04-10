# CMake >= 3.28 suggested
# CMake >= 3.25 required for fine-grained scope with block() and enabling configure time tests
# CMake >= 3.24: find_package(... GLOBAL) and FetchContent_Declare(... FIND_PACKAGE_ARGS)
# CMake >= 3.21:  *_IS_TOP_LEVEL variable used to manage subproject dependencies
cmake_minimum_required(VERSION 3.25...4.3)

# --- CMAKE_BUILD_TYPE default
# The simulations are 10x slower for default to Debug.
# Thus, for single config generators, set build type to Release
get_property(is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT is_multi_config AND NOT (CMAKE_BUILD_TYPE OR DEFINED ENV{CMAKE_BUILD_TYPE}))
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Release can be 10x faster simulation run time for gemini3d.run vs. Debug")
endif()

if(CMAKE_GENERATOR MATCHES "^Visual Studio")
  message(FATAL_ERROR "${CMAKE_GENERATOR} is not supported. Use \"cmake -G Ninja\" option.")
endif()

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.26.0 AND CMAKE_VERSION VERSION_LESS 3.28.0)
  message(WARNING "CMake 3.26.* and 3.27.* have known bugs causing build failures due to incorrectly graphing Fortran module dependencies, regardless of operating system or compiler.
Suggest CMake 3.28 or newer https://github.com/Kitware/CMake/releases
If on Linux HPC, do like 'module avail cmake' to see if a suitable CMake version is available.")
endif()

# --- main Gemini3D build

project(gemini3d
LANGUAGES C CXX Fortran
# Gemini3D is Fortran, but external libraries use C, and some find_package need C.
DESCRIPTION "3-D ionospheric model"
HOMEPAGE_URL https://github.com/gemini3d/gemini
VERSION 2.0.0
)

enable_testing()  # keep this so BUILD_TESTING=off doesn't remove all tests
include(CTest)

include(FetchContent)

set(CMAKE_CXX_SCAN_FOR_MODULES OFF)
# This project isn't using C++ modules; avoid bugs in module scanning CMake 3.28, 3.29 with GCC >= 14.

include(cmake/package/git_rev.cmake)
include(options.cmake)

set(CMAKE_CXX_STANDARD 17)

file(READ cmake/libraries.json json)

# Fortran filesystem library "ffilesystem"
set(ffilesystem_fallback true)
set(ffilesystem_extra false)
set(HAVE_F03TYPE false)
# may need fallback for old or slightly broken compilers
string(JSON ffilesystem_url GET "${json}" "ffilesystem")
FetchContent_Declare(ffilesystem URL ${ffilesystem_url})
FetchContent_MakeAvailable(ffilesystem)

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

string(JSON h5fortran_url GET "${json}" "h5fortran")
FetchContent_Declare(h5fortran URL ${h5fortran_url})
FetchContent_MakeAvailable(h5fortran)
# this also implicitly finds HDF5.
# h5fortran is a high-level, object-oriented HDF5 interface.

include(cmake/compilers.cmake)

# --- linear algebra libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
if(NOT DEFINED ENV{MKLROOT} AND NOT LAPACK_VENDOR MATCHES "^MKL")
  # oneMKL MKLConfig.cmake must be invoked only once, and MUMPS invokes it.
  find_package(LAPACK REQUIRED GLOBAL)
else()
  find_package(SCALAPACK REQUIRED GLOBAL)
endif()

set(BUILD_SINGLE off)
set(BUILD_DOUBLE on)
set(BUILD_COMPLEX off)
set(BUILD_COMPLEX16 off)

string(JSON mumps_url GET "${json}" "mumps")
FetchContent_Declare(mumps URL ${mumps_url})
FetchContent_MakeAvailable(mumps)

# --- climate models
if(gemini3d_glow)
  string(JSON glow_url GET "${json}" "glow")
  FetchContent_Declare(glow URL ${glow_url})
  FetchContent_MakeAvailable(glow)
endif()

if(gemini3d_hwm14)
  string(JSON hwm14_url GET "${json}" "hwm14")
  FetchContent_Declare(hwm14 URL ${hwm14_url})
  FetchContent_MakeAvailable(hwm14)
endif()

# --- MSISE00 / MSIS 2.x
set(msis_BUILD_UTILS on)
string(JSON msis_url GET "${json}" "msis")
FetchContent_Declare(msis URL ${msis_url})
FetchContent_MakeAvailable(msis)

# --- Gemini3D library
add_subdirectory(src)

# --- Gemini3D self test
include(${CMAKE_CURRENT_SOURCE_DIR}/test/mpi_launcher.cmake) # for MPI tests

add_subdirectory(test)

# summary print
include(cmake/summary.cmake)

# packaging
include(cmake/package/pkgconf.cmake)
include(cmake/package/install.cmake)
