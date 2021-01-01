# Finds Lapack, tests, and if not found or broken, autobuild Lapack

if(autobuild)
  find_package(LAPACK)
else()
  find_package(LAPACK REQUIRED)
endif()

if(NOT LAPACK_FOUND)
  set(lapack_external true CACHE BOOL "autobuild Lapack")

  include(FetchContent)

  FetchContent_Declare(lapack_proj
    GIT_REPOSITORY ${gemini_lapack_url}
    GIT_TAG ${gemini_lapack_tag}
    GIT_SHALLOW true
    CMAKE_ARGS "-Darith=${arith}"
  )

  FetchContent_MakeAvailable(lapack_proj)

  add_library(LAPACK::LAPACK ALIAS lapack)
  add_library(BLAS::BLAS ALIAS blas)
endif()
