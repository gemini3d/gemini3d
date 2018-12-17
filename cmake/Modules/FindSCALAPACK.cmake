# https://github.com/certik/hermes/blob/master/hermes_common/cmake/FindSCALAPACK.cmake
# ScaLAPACK and BLACS

# USEMKL: Using MKL with GNU or other compiler

unset(SCALAPACK_LIBRARY)

if(USEMKL OR CMAKE_Fortran_COMPILER_ID STREQUAL Intel)

  # FIXME: this would be for threaded
  # mkl_scalapack_lp64 mkl_intel_lp64 mkl_intel_thread mkl_core mkl_blacs_intelmpi_lp64 iomp5
  
  if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
    set(slibs mkl_scalapack_lp64 mkl_intel_lp64 mkl_sequential mkl_core mkl_blacs_intelmpi_lp64)
  else()
    set(slibs mkl_scalapack_lp64 mkl_gf_lp64 mkl_sequential mkl_core mkl_blacs_openmpi_lp64)
  endif()

  # this is for sequential:
  foreach(slib ${slibs})
    find_library(SCALAPACK_${slib}_LIBRARY
             NAMES ${slib}
             PATHS $ENV{MKLROOT}/lib
                   $ENV{MKLROOT}/lib/intel64
                   $ENV{INTEL}/mkl/lib/intel64
             NO_DEFAULT_PATH)
    if(NOT SCALAPACK_${slib}_LIBRARY)
      message(FATAL_ERROR "NOT FOUND: " ${slib})
    endif()
    message(STATUS "Intel MKL Scalapack FOUND: " ${slib})
    list(APPEND SCALAPACK_LIBRARY ${SCALAPACK_${slib}_LIBRARY})
    mark_as_advanced(SCALAPACK_${slib}_LIBRARY)
  endforeach()
  list(APPEND SCALAPACK_LIBRARY pthread ${CMAKE_DL_LIBS} m)

else()

    find_library(SCALAPACK_LIBRARY
                 NAMES scalapack scalapack-pvm scalapack-mpi scalapack-mpich scalapack-mpich2 scalapack-openmpi scalapack-lam
                 PATH_SUFFIXES lib)

endif()


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SCALAPACK
    REQUIRED_VARS SCALAPACK_LIBRARY)  
# don't put BLACS_LIBRARY REQUIRED_VARS because it might be in libscalapack.a)

if(SCALAPACK_FOUND)
  set(SCALAPACK_LIBRARIES ${SCALAPACK_LIBRARY})
endif()

mark_as_advanced(SCALAPACK_LIBRARY)
