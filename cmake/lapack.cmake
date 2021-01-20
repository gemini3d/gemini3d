# Finds Lapack, tests, and if not found or broken, autobuild Lapack

if(autobuild)
  find_package(LAPACK)
else()
  find_package(LAPACK REQUIRED)
endif()

if(NOT LAPACK_FOUND)
  set(lapack_external true CACHE BOOL "autobuild Lapack")

  include(FetchContent)

  if(GIT_FOUND)
    FetchContent_Declare(lapack_proj
      GIT_REPOSITORY ${lapack_git}
      GIT_TAG ${lapack_tag}
      GIT_SHALLOW true
      CMAKE_ARGS "-Darith=${arith}"
    )
  else(GIT_FOUND)
    FetchContent_Declare(lapack_proj
      URL ${lapack_zip}
      TLS_VERIFY true
      UPDATE_DISCONNECTED true
      CMAKE_ARGS "-Darith=${arith}"
    )
  endif(GIT_FOUND)

  FetchContent_MakeAvailable(lapack_proj)

  add_library(LAPACK::LAPACK ALIAS lapack)
  add_library(BLAS::BLAS ALIAS blas)
endif()
