cmake_minimum_required(VERSION 3.17...3.25)
# CMake presets require >= 3.21
# full test capabilities required >= 3.22

# --- CMAKE_BUILD_TYPE default
# Normally we shouldn't set a default like this, but the simulations are 10x slower
# if the user lets it default to Debug.
# for single config generators, set build type to Release
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

message(STATUS "Gemini3D: CMake ${CMAKE_VERSION}")

if(NOT DEFINED ${PROJECT_NAME}_BUILD_TESTING)
  set(${PROJECT_NAME}_BUILD_TESTING ${BUILD_TESTING})
endif()

include(cmake/package/git_rev.cmake)
include(cmake/options.cmake)

# check that gemini3d/external libraries are installed
set(need_gemext "
Gemini3D requires several external libraries that are one-time installed via this procedure.
'~/libgem' directory is an arbitrary location.

  git clone https://github.com/gemini3d/external
  cmake -S external -B external/build -DCMAKE_INSTALL_PREFIX=~/libgem
  cmake --build external/build

Now, configure and build Gemini3D (from gemini3d/ directory) like this:

  cmake -B build -DCMAKE_PREFIX_PATH=~/libgem
  cmake --build build --parallel
")

if(NOT IS_DIRECTORY ${CMAKE_PREFIX_PATH}/cmake)
  message(FATAL_ERROR ${need_gemext})
endif()

# Fortran filesystem library "ffilesystem"
find_package(ffilesystem CONFIG)
# this is the first package to find, so make a special error message for users
if(NOT ffilesystem_FOUND)
  message(FATAL_ERROR ${need_gemext})
endif()


find_package(MPI COMPONENTS C CXX Fortran REQUIRED)

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
if(CMAKE_VERSION VERSION_LESS 3.20)
  get_filename_component(msis_dir "${msis_loc}" DIRECTORY)
else()
  cmake_path(GET msis_loc PARENT_PATH msis_dir)
endif()
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
