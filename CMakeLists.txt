cmake_minimum_required(VERSION 3.15...3.26)

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
include(cmake/options.cmake)

# Fortran filesystem library "ffilesystem"
find_package(ffilesystem CONFIG)
# this is the first package to find, so make a special error message for users
if(NOT ffilesystem_FOUND)
  message(FATAL_ERROR ${need_gemext})
endif()
if(ffilesystem_VERSION VERSION_LESS 4.0)
  message(STATUS "Skipping C++ frontend due to ffilesystem version ${ffilesystem_VERSION} < 4.0")
  set(cpp false)
endif()

# for gemini_c C++ front end
if(cpp)
  find_package(iniparser CONFIG)
  if(NOT iniparser_FOUND)
    set(cpp false)
  endif()
endif()

# --- MPI
# NOTE: our find_package(MPI) needs to be before find_package(MUMPS), which also calls find_package(MPI)
# Conda e.g. scikit-learn can cause problems with finding MPI, so exclude Conda from search
if(DEFINED ENV{CONDA_PREFIX})
  set(mpi_ignore_path
    $ENV{CONDA_PREFIX}/bin $ENV{CONDA_PREFIX}/lib $ENV{CONDA_PREFIX}/include
    $ENV{CONDA_PREFIX}/Library/bin $ENV{CONDA_PREFIX}/Library/lib $ENV{CONDA_PREFIX}/Library/include
  )
  list(APPEND CMAKE_IGNORE_PATH ${mpi_ignore_path})
endif()
set(MPI_DETERMINE_LIBRARY_VERSION true)
find_package(MPI COMPONENTS C CXX Fortran REQUIRED)
message(STATUS "MPI Library Version: ${MPI_C_LIBRARY_VERSION_STRING}")
include(cmake/openmpi.cmake)

# --- end MPI

find_package(h5fortran CONFIG REQUIRED)

include(cmake/compilers.cmake)

# --- linear algebra libraries
find_package(LAPACK REQUIRED)

find_package(MUMPS REQUIRED COMPONENTS d)

# --- climate models
if(glow)
  find_package(glow CONFIG REQUIRED)
endif()

if(hwm14)
  find_package(hwm14 CONFIG REQUIRED)
endif()

# --- MSISE00 or MSIS 2.0
find_package(msis CONFIG REQUIRED)
# WORKAROUND: MSIS shared Rpath fail
# Rpath is defined in MSIS, but ld fails to use it
get_target_property(msis_loc msis::msis_ifc LOCATION)
get_filename_component(msis_dir "${msis_loc}" DIRECTORY)

# END WORKAROUND

# pop off ignored paths so rest of script can find Python
if(DEFINED CMAKE_IGNORE_PATH)
  list(REMOVE_ITEM CMAKE_IGNORE_PATH "${mpi_ignore_path}")
endif()


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
