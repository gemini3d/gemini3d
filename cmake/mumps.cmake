# MUMPS -- use cmake -DMUMPS_ROOT= for hint.
#
# Intel MKL-compiled MUMPS requires at the linker for the main executable:
# mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 mkl_intel_lp64 mkl_intel_thread mkl_core
#
# easily obtain MUMPS without compiling:
# CentOS 6/7 EPEL: yum install mumps-devel
# Ubuntu / Debian: apt install libmumps-dev

unset(_mumps_extra)

include(${CMAKE_CURRENT_LIST_DIR}/scalapack.cmake)

# --- MUMPS
if(realbits EQUAL 32)
  set(arith s)
else()
  set(arith d)
endif()

set(mumps_external true)
if(MUMPS_ROOT OR CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  find_package(MUMPS COMPONENTS ${arith})
  if(MUMPS_FOUND)
    set(mumps_external false)
  endif()
endif()
if(mumps_external)
  include(${CMAKE_CURRENT_LIST_DIR}/mumps_external.cmake)
endif()

if(metis)
  find_package(METIS REQUIRED)
  list(APPEND _mumps_extra ${METIS_LIBRARIES})
endif()
if(scotch)
  find_package(Scotch REQUIRED COMPONENTS ESMUMPS)
  list(APPEND _mumps_extra ${Scotch_LIBRARIES})
endif()
# rather than appending libraries everywhere, just put them together here.
list(APPEND MUMPS_LIBRARIES ${SCALAPACK_LIBRARIES} ${LAPACK_LIBRARIES} ${_mumps_extra})
if(OpenMP_FOUND)
  list(APPEND MUMPS_LIBRARIES OpenMP::OpenMP_Fortran OpenMP::OpenMP_C)
endif()
list(APPEND MUMPS_INCLUDE_DIRS ${SCALAPACK_INCLUDE_DIRS})

if(mumps_external)
  return()
endif()

# -- minimal check that MUMPS is linkable
set(CMAKE_REQUIRED_LIBRARIES ${MUMPS_LIBRARIES} MPI::MPI_Fortran)
set(CMAKE_REQUIRED_INCLUDES ${MUMPS_INCLUDE_DIRS})

check_fortran_source_compiles("include '${arith}mumps_struc.h'
type(${arith}mumps_struc) :: mumps_par
end"
  MUMPS_OK SRC_EXT f90)

if(NOT MUMPS_OK)
message(FATAL_ERROR "MUMPS ${MUMPS_LIBRARIES} not working with ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
endif()
