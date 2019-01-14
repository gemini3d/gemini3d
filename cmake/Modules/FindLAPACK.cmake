# LAPACK and BLAS
# CMake factory FindLAPACK at least through 3.13 ignores LAPACK_ROOT and is a mess, so we made this overly simple FindLAPACK.

cmake_policy(VERSION 3.3)

find_library(LAPACK_LIBRARY
  NAMES lapack)

find_library(BLAS_LIBRARY
  NAMES blas)            


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  LAPACK
  REQUIRED_VARS LAPACK_LIBRARY BLAS_LIBRARY)

if(SCALAPACK_FOUND)
  set(LAPACK_LIBRARIES ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
endif()

mark_as_advanced(LAPACK_LIBRARY BLAS_LIBRARY)

