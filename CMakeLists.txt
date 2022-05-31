cmake_minimum_required(VERSION 3.20...3.23)
# CMake presets require >= 3.21

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
include(CTest)

message(STATUS "Gemini3D: CMake ${CMAKE_VERSION}")

if(NOT DEFINED ${PROJECT_NAME}_BUILD_TESTING)
  set(${PROJECT_NAME}_BUILD_TESTING ${BUILD_TESTING})
endif()

include(cmake/package/git_rev.cmake)
include(cmake/config/options.cmake)

# Fortran filesystem library "ffilesystem"
find_package(ffilesystem CONFIG)
# this is the first package to find, so make a special error message for users
if(NOT ffilesystem_FOUND OR ffilesystem_VERSION VERSION_LESS 3.3.0)
  message(FATAL_ERROR "
Gemini3D requires several external libraries that are one-time installed via this procedure. 'gemini_libs' directory is an arbitrary location.

  git clone https://github.com/gemini3d/external
  cmake -S external -B external/build -DCMAKE_INSTALL_PREFIX=~/gemini_libs
  cmake --build external/build

Now, configure and build Gemini3D (from gemini3d/ directory) like this:

  cmake -B build -DCMAKE_PREFIX_PATH=~/gemini_libs
  cmake --build build --parallel
")
endif()
if(ffilesystem_VERSION VERSION_LESS 3.4.2 AND CMAKE_CXX_COMPILER_ID STREQUAL GNU AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.1.0)
  #target_link_libraries(ffilesystem::filesystem INTERFACE stdc++fs)
  # function target_link_libraries() should work, but doesn't (CMake bug)
  set_target_properties(ffilesystem::filesystem PROPERTIES INTERFACE_LINK_LIBRARIES stdc++fs)
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

# --- MSISE00 or MSIS 2.0
find_package(msis CONFIG REQUIRED)
# WORKAROUND: MSIS shared Rpath fail
# Rpath is defined in MSIS, but ld fails to use it
get_target_property(msis_loc msis::msis_ifc LOCATION)
cmake_path(GET msis_loc PARENT_PATH msis_dir)
# END WORKAROUND


include(cmake/compilers/flags.cmake)

# unit test DLL paths
include(cmake/DllTestPath.cmake)


add_subdirectory(src)
# Gemini3D source code

add_subdirectory(test)
# fundamental tests of MPI and numeric libraries essential for Gemini3D

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
