# Finds Scalapack, tests, and if not found or broken, autobuild scalapack
include(FetchContent)

if(NOT scalapack_external)
  if(autobuild)
    find_package(SCALAPACK)
  else()
    find_package(SCALAPACK REQUIRED)
  endif()
endif()

if(SCALAPACK_FOUND OR TARGET SCALAPACK::SCALAPACK)
  return()
endif()

set(scalapack_external true CACHE BOOL "build ScaLapack")

FetchContent_Declare(SCALAPACK
  GIT_REPOSITORY ${scalapack_git}
  GIT_TAG ${scalapack_tag}
  CMAKE_ARGS -Darith=${arith})

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.14)
  FetchContent_MakeAvailable(SCALAPACK)
elseif(NOT scalapack_POPULATED)
  FetchContent_Populate(SCALAPACK)
  add_subdirectory(${scalapack_SOURCE_DIR} ${scalapack_BINARY_DIR})
endif()
