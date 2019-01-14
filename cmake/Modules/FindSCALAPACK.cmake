# https://github.com/certik/hermes/blob/master/hermes_common/cmake/FindSCALAPACK.cmake
# ScaLAPACK and BLACS

# USEMKL: Using MKL with GNU or other compiler

cmake_policy(VERSION 3.3)

unset(SCALAPACK_LIBRARY)
unset(SCALAPACK_OpenMPI_FOUND)
unset(SCALAPACK_MPICH_FOUND)

if(NOT SCALAPACK_FIND_COMPONENTS)
  set(SCALAPACK_FIND_COMPONENTS OpenMPI)
endif()

function(mkl_scala)

foreach(s ${ARGV})
  find_library(SCALAPACK_${s}_LIBRARY
           NAMES ${s}
           PATHS $ENV{MKLROOT}/lib
                 $ENV{MKLROOT}/lib/intel64
                 $ENV{INTEL}/mkl/lib/intel64
           NO_DEFAULT_PATH)
  if(NOT SCALAPACK_${s}_LIBRARY)
    message(FATAL_ERROR "NOT FOUND: " ${s})
  endif()
  
  list(APPEND SCALAPACK_LIB ${SCALAPACK_${s}_LIBRARY})
endforeach()

list(APPEND SCALAPACK_LIB pthread ${CMAKE_DL_LIBS} m)

set(SCALAPACK_LIBRARY ${SCALAPACK_LIB} PARENT_SCOPE)

endfunction()



if(IntelPar IN_LIST SCALAPACK_FIND_COMPONENTS)

  mkl_scala(mkl_scalapack_lp64 mkl_intel_lp64 mkl_intel_thread mkl_core mkl_blacs_intelmpi_lp64 iomp5)

  if(SCALAPACK_LIBRARY)
    set(SCALAPACK_IntelPar_FOUND true)
  endif()

elseif(IntelSeq IN_LIST SCALAPACK_FIND_COMPONENTS)

  mkl_scala(mkl_scalapack_lp64 mkl_intel_lp64 mkl_sequential mkl_core mkl_blacs_intelmpi_lp64)

  if(SCALAPACK_LIBRARY)
    set(SCALAPACK_IntelSeq_FOUND true)
  endif()

elseif(MKL IN_LIST SCALAPACK_FIND_COMPONENTS)

  mkl_scala(mkl_scalapack_lp64 mkl_gf_lp64 mkl_sequential mkl_core mkl_blacs_openmpi_lp64)

  if(SCALAPACK_LIBRARY)
    set(SCALAPACK_MKL_FOUND true)
  endif()
  
elseif(OpenMPI IN_LIST SCALAPACK_FIND_COMPONENTS)

  find_library(SCALAPACK_LIBRARY
               NAMES scalapack scalapack-openmpi
               PATH_SUFFIXES lib)
               
  if(SCALAPACK_LIBRARY)
    set(SCALAPACK_OpenMPI_FOUND true)
  endif()

elseif(MPICH IN_LIST SCALAPACK_FIND_COMPONENTS)
  find_library(SCALAPACK_LIBRARY
               NAMES scalapack-mpich scalapack-mpich2
               PATH_SUFFIXES lib)
               
  if(SCALAPACK_LIBRARY)
    set(SCALAPACK_MPICH_FOUND true)
  endif()

endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  SCALAPACK
  REQUIRED_VARS SCALAPACK_LIBRARY
  HANDLE_COMPONENTS)

if(SCALAPACK_FOUND)
  set(SCALAPACK_LIBRARIES ${SCALAPACK_LIBRARY})
endif()

mark_as_advanced(SCALAPACK_LIBRARY)

