set(lapack_external true CACHE BOOL "autobuild Lapack")

include(FetchContent)

FetchContent_Declare(lapack_proj
  GIT_REPOSITORY https://github.com/scivision/lapack.git
  GIT_TAG v3.9.0.2
  CMAKE_ARGS "-Darith=${arith}"
)

FetchContent_MakeAvailable(lapack_proj)

add_library(LAPACK::LAPACK ALIAS lapack)
add_library(BLAS::BLAS ALIAS blas)
set(LAPACK_FOUND true)
set(BLAS_FOUND true)
