# MUMPS -- use cmake -DMUMPS_ROOT= for hint.
#
# Intel MKL-compiled MUMPS requires at the linker for the main executable:
# mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 mkl_intel_lp64 mkl_intel_thread mkl_core
#
# easily obtain MUMPS without compiling:
# CentOS 6/7 EPEL: yum install mumps-devel
# Ubuntu / Debian: apt install libmumps-dev

# --- prereqs
include(${CMAKE_CURRENT_LIST_DIR}/lapack.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/scalapack.cmake)

# --- MUMPS

if(mumps_external)
  include(${CMAKE_CURRENT_LIST_DIR}/mumps_external.cmake)
  return()
endif()

unset(_mumps_extra)

if(MUMPS_ROOT OR (DEFINED ENV{MUMPS_ROOT}) OR (CMAKE_Fortran_COMPILER_ID STREQUAL GNU))
  find_package(MUMPS COMPONENTS ${arith})
else()
  message(VERBOSE "Skipping find_package(MUMPS)")
endif()
if(NOT MUMPS_FOUND)
  include(${CMAKE_CURRENT_LIST_DIR}/mumps_external.cmake)
  set(mumps_external true CACHE BOOL "autobuild Mumps")
else()
  set(mumps_external false CACHE BOOL "autobuild Mumps")
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
list(APPEND MUMPS_INCLUDE_DIRS ${SCALAPACK_INCLUDE_DIRS} ${LAPACK_INCLUDE_DIRS})

if(mumps_external OR scalapack_external OR lapack_external)
# pre-build checks can't be used when external library isn't built yet.
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
  message(STATUS "MUMPS ${MUMPS_LIBRARIES} not working with ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
  include(${CMAKE_CURRENT_LIST_DIR}/mumps_external.cmake)
  set(mumps_external true CACHE BOOL "autobuild Mumps")
endif()
