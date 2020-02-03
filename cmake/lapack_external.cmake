include(FetchContent)

FetchContent_Declare(lapack_proj
  GIT_REPOSITORY https://github.com/scivision/lapack.git
  GIT_TAG v3.9.0.1
  CMAKE_ARGS "-Darith=${arith}"
)

FetchContent_MakeAvailable(lapack_proj)

if(NOT LAPACK_FOUND)
  message(STATUS "AUTOBUILD: LAPACK + BLAS")
endif()

set(LAPACK_LIBRARIES lapack)
set(BLAS_LIBRARIES blas)
set(LAPACK_FOUND true)
set(BLAS_FOUND true)