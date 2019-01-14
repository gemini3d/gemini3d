# LAPACK and BLAS
# CMake factory FindLAPACK at least through 3.13 ignores LAPACK_ROOT and is a mess, so we made this overly simple FindLAPACK.

cmake_policy(VERSION 3.3)

function(mkl_libs)
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor

foreach(s ${ARGV})
  find_library(LAPACK_${s}_LIBRARY
           NAMES ${s}
           PATHS $ENV{MKLROOT}/lib
                 $ENV{MKLROOT}/lib/intel64
                 $ENV{MKLROOT}/../compiler/lib/intel64
           NO_DEFAULT_PATH)
  if(NOT LAPACK_${s}_LIBRARY)
    message(FATAL_ERROR "NOT FOUND: " ${s})
  endif()
  
  list(APPEND LAPACK_LIB ${LAPACK_${s}_LIBRARY})
endforeach()

list(APPEND LAPACK_LIB pthread ${CMAKE_DL_LIBS} m)

set(LAPACK_LIBRARY ${LAPACK_LIB} PARENT_SCOPE)
set(BLAS_LIBRARY ${LAPACK_LIB} PARENT_SCOPE)

endfunction()

if(USEMKL OR CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  mkl_libs(mkl_intel_lp64 mkl_intel_thread mkl_core iomp5)
else()
  find_library(LAPACK_LIBRARY
    NAMES lapack)

  find_library(BLAS_LIBRARY
    NAMES blas)            
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  LAPACK
  REQUIRED_VARS LAPACK_LIBRARY BLAS_LIBRARY)

if(SCALAPACK_FOUND)
  set(LAPACK_LIBRARIES ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
endif()

mark_as_advanced(LAPACK_LIBRARY BLAS_LIBRARY)

