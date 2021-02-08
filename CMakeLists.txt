cmake_minimum_required(VERSION 3.13)

include(cmake/policy.cmake)
include(cmake/compiler_find.cmake)
find_c_fortran()

file(STRINGS cmake/gemini_version.txt PROJECT_VERSION REGEX "^([0-9]+\.[0-9]+\.[0-9]+)" LIMIT_INPUT 16 LENGTH_MAXIMUM 16 LIMIT_COUNT 1)

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

include(cmake/git_rev.cmake)
include(cmake/options.cmake)
include(cmake/libraries.cmake)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)

include(cmake/h5fortran.cmake)
include(cmake/nc4fortran.cmake)

# do these after h5fortran to avoid threads test failure
include(cmake/compilers.cmake)
# first this to get libs
include(cmake/mumps.cmake)
include(cmake/compiler_flags.cmake)
# this comes last to avoid intermittent failures

# all Gemini3D source code under this
add_subdirectory(src)

# self-test simulations -- after all targets for if(TARGET ...)
include(cmake/matlab.cmake)
include(cmake/python.cmake)
include(cmake/test_setup.cmake)
include(cmake/test_sim.cmake)

# summary print
include(cmake/summary.cmake)
