# Finds Scalapack, tests, and if not found or broken, autobuild scalapack
include(FetchContent)

if(autobuild)
  find_package(SCALAPACK)
else()
  find_package(SCALAPACK REQUIRED)
endif()

if(NOT SCALAPACK_FOUND)
  set(scalapack_external true CACHE BOOL "autobuild ScaLapack")

  FetchContent_Declare(SCALAPACK
    GIT_REPOSITORY ${scalapack_git}
    GIT_TAG ${scalapack_tag}
    CMAKE_ARGS "-Darith=${arith}")

  FetchContent_MakeAvailable(SCALAPACK)
endif()

target_link_libraries(SCALAPACK::SCALAPACK INTERFACE LAPACK::LAPACK)
