cmake_minimum_required(VERSION 3.19...3.21)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Debug or Release")
endif()

file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/VERSION PROJECT_VERSION
  REGEX "^([0-9]+\.[0-9]+\.[0-9]+)" LIMIT_INPUT 16 LENGTH_MAXIMUM 16 LIMIT_COUNT 1)

project(gemini3d
  LANGUAGES C Fortran
  # Gemini3D is Fortran, but external libraries use C, and some find_package need C.
  DESCRIPTION "3-D ionospheric model"
  HOMEPAGE_URL https://github.com/gemini3d/gemini
  VERSION ${PROJECT_VERSION})

enable_testing()  # keep this so BUILD_TESTING=off doesn't remove all tests
include(CTest) # for CDash
if(NOT DEFINED ${PROJECT_NAME}_BUILD_TESTING)
  set(${PROJECT_NAME}_BUILD_TESTING ${BUILD_TESTING})
endif()

include(cmake/package/git_rev.cmake)
include(cmake/config/options.cmake)
include(cmake/config/libraries.cmake)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)

include(cmake/ext_libs/h5fortran.cmake)
include(cmake/ext_libs/nc4fortran.cmake)

# do these after h5fortran to avoid threads test failure
include(cmake/compilers/compilers.cmake)
# first this to get libs
include(cmake/ext_libs/mumps.cmake)
include(cmake/compilers/flags.cmake)
# flags.cmake comes last to avoid intermittent failures

# utility
include(cmake/SetTargProps.cmake)
# unit test DLL paths
include(cmake/test/lib_path.cmake)

# all Gemini3D source code under this
add_subdirectory(src)

if(CMAKE_VERSION VERSION_LESS 3.20)
  message(STATUS "*** Gemini3D self-tests require CMake >= 3.20 ***
  Recommendation: Install latest CMake by:
    cmake -P ${PROJECT_SOURCE_DIR}/scripts/install_cmake.cmake")
  return()
endif()

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
