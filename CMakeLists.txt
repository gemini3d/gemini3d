cmake_minimum_required(VERSION 3.15)
# 3.14 for check_fortran_source_runs & fetchcontent
# 3.15 for robust Python finding

include(cmake/policy.cmake)

project(gemini3d
  LANGUAGES C Fortran  # MUST include C language for Intel / MKL to work
  DESCRIPTION "3-D ionospheric model"
  HOMEPAGE_URL https://github.com/gemini3d/gemini
  VERSION 0.6.1)

enable_testing()  # keep this so BUILD_TESTING=off doesn't remove all tests
include(CTest) # for CDash
if(NOT DEFINED ${PROJECT_NAME}_BUILD_TESTING)
  set(${PROJECT_NAME}_BUILD_TESTING ${BUILD_TESTING})
endif()

include(cmake/options.cmake)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/h5fortran.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/nc4fortran.cmake)

# do these after h5fortran to avoid h5fortran HDF5OK threads test failure
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/compilers.cmake)
# first this to get libs
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/mumps.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/compiler_flags.cmake)
# this comes last to avoid intermittent failures

# self-test simulations
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/python.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/test_setup.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/test_sim.cmake)

add_subdirectory(src)

# summary print
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/summary.cmake)
