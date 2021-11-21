cmake_minimum_required(VERSION 3.20...3.22)

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
LANGUAGES C Fortran
# Gemini3D is Fortran, but external libraries use C, and some find_package need C.
DESCRIPTION "3-D ionospheric model"
HOMEPAGE_URL https://github.com/gemini3d/gemini
VERSION ${PROJECT_VERSION}
)

enable_testing()  # keep this so BUILD_TESTING=off doesn't remove all tests
include(CTest) # for CDash

if(NOT DEFINED ${PROJECT_NAME}_BUILD_TESTING)
  set(${PROJECT_NAME}_BUILD_TESTING ${BUILD_TESTING})
endif()

if(CMAKE_GENERATOR STREQUAL "Ninja Multi-Config")
  set(EXTPROJ_GENERATOR "Ninja")
else()
  set(EXTPROJ_GENERATOR ${CMAKE_GENERATOR})
endif()

include(cmake/package/git_rev.cmake)
include(cmake/config/options.cmake)
include(cmake/config/libraries.cmake)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)

if(mpi)
  find_package(MPI COMPONENTS C Fortran REQUIRED)
else()
  add_subdirectory(src/vendor/mpi_stubs)
endif()

include(cmake/ext_libs/h5fortran.cmake)
include(cmake/ext_libs/nc4fortran.cmake)

# do these after h5fortran to avoid threads test failure
include(cmake/compilers/compilers.cmake)
if(mpi)
  # here for HWLOC_FOUND directory scope
  find_package(HWLOC)
endif()
# first this to get libs
include(cmake/ext_libs/mumps.cmake)
if(glow)
  include(${PROJECT_SOURCE_DIR}/cmake/ext_libs/glow.cmake)
endif(glow)

include(cmake/compilers/flags.cmake)
# flags.cmake comes last to avoid intermittent failures

# unit test DLL paths
include(cmake/test/lib_path.cmake)

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
