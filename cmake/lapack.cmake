# Finds Lapack, tests, and if not found or broken, autobuild Lapack
include(FetchContent)

if(NOT lapack_external)
  if(autobuild)
    find_package(LAPACK)
  else()
    find_package(LAPACK REQUIRED)
  endif()
endif()

if(LAPACK_FOUND)
  return()
endif()


set(lapack_external true CACHE BOOL "build Lapack")

FetchContent_Declare(LAPACK
  GIT_REPOSITORY ${lapack_git}
  GIT_TAG ${lapack_tag}
  CMAKE_ARGS "-Darith=${arith}")

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.14)
  FetchContent_MakeAvailable(LAPACK)
elseif(NOT lapack_POPULATED)
  FetchContent_Populate(LAPACK)
  add_subdirectory(${lapack_SOURCE_DIR} ${lapack_BINARY_DIR})
endif()


add_library(LAPACK::LAPACK ALIAS lapack)
add_library(BLAS::BLAS ALIAS blas)
