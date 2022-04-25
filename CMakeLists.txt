cmake_minimum_required(VERSION 3.20...3.23)
# For CMake presets, CMake >= 3.21 required

file(READ ${CMAKE_CURRENT_LIST_DIR}/codemeta.json _j)
string(JSON PROJECT_VERSION GET ${_j} version)

# --- CMAKE_BUILD_TYPE default
# Normally we shouldn't set a default like this, but the simulations are 10x slower
# if the user lets it default to Debug.
# for single config generators, set build type to Release
get_property(is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT is_multi_config AND NOT (CMAKE_BUILD_TYPE OR DEFINED ENV{CMAKE_BUILD_TYPE}))
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Release can be 10x faster simulation run time for gemini3d.run vs. Debug")
endif()

# --- main Gemini3D build

project(gemini3d
LANGUAGES C CXX Fortran
# Gemini3D is Fortran, but external libraries use C, and some find_package need C.
DESCRIPTION "3-D ionospheric model"
HOMEPAGE_URL https://github.com/gemini3d/gemini
VERSION ${PROJECT_VERSION}
)

enable_testing()  # keep this so BUILD_TESTING=off doesn't remove all tests
include(CTest) # for CDash

message(STATUS "Gemini3D: CMake ${CMAKE_VERSION}")

if(NOT DEFINED ${PROJECT_NAME}_BUILD_TESTING)
  set(${PROJECT_NAME}_BUILD_TESTING ${BUILD_TESTING})
endif()

include(cmake/package/git_rev.cmake)
include(cmake/config/options.cmake)

# Fortran filesystem library "ffilesystem"
find_package(ffilesystem CONFIG)
if(ffilesystem_VERSION VERSION_LESS 3.4.2)
  target_link_libraries(ffilesystem::filesystem INTERFACE
  $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS:${CMAKE_CXX_COMPILER_VERSION},9.1.0>>:stdc++fs>
  )
endif()
# this is the first package we find, so make a special error message for users
if(NOT ffilesystem_FOUND OR ffilesystem_VERSION VERSION_LESS 3.3.0)
  message(FATAL_ERROR "
Gemini3D requires several external libraries that are one-time installed via this procedure. Note: the 'gemini_libs' directory is arbitrary example.

  git clone https://github.com/gemini3d/external
  cmake -B external/build -S external -DCMAKE_INSTALL_PREFIX=~/gemini_libs
  cmake --build external/build

Now, configure and build Gemini3D (from gemini3d/ directory) like this:

  cmake -B build -DCMAKE_PREFIX_PATH=~/gemini_libs
  cmake --build build --parallel
")
endif()

if(mpi)
  find_package(MPI COMPONENTS C CXX Fortran REQUIRED)
else()
  find_package(mpi_scalapack_stub CONFIG REQUIRED)
endif()

if(hdf5)
  find_package(h5fortran CONFIG REQUIRED)
endif()
if(netcdf)
  find_package(nc4fortran CONFIG REQUIRED)
endif()

include(cmake/compilers/compilers.cmake)
if(mpi)
  # here for HWLOC_FOUND directory scope
  find_package(HWLOC)
endif()

# --- linear algebra libraries
find_package(LAPACK REQUIRED)
if(mpi)
  find_package(MUMPS REQUIRED COMPONENTS d)
else()
  find_package(MUMPS REQUIRED COMPONENTS d mpiseq)
endif()

# --- climate models
if(glow)
  find_package(glow CONFIG REQUIRED)
endif()

if(hwm14)
  find_package(hwm14 CONFIG REQUIRED)
endif()

include(cmake/compilers/flags.cmake)

# unit test DLL paths
include(cmake/DllTestPath.cmake)

# all Gemini3D source code under this
add_subdirectory(src)

# self-test simulations -- after all targets for if(TARGET ...)
if(matlab)
  include(cmake/ext_libs/matlab.cmake)
endif(matlab)
if(python)
  include(cmake/ext_libs/python.cmake)
endif(python)

include(cmake/test/config.cmake)
include(cmake/test/sim.cmake)

# summary print
include(cmake/config/summary.cmake)

# packaging
include(cmake/package/pkgconf.cmake)
include(cmake/package/install.cmake)
