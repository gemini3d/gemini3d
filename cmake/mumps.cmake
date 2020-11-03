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

if(mpi)
  include(${CMAKE_CURRENT_LIST_DIR}/scalapack.cmake)
endif()
# --- MUMPS

if(MUMPS_ROOT OR (DEFINED ENV{MUMPS_ROOT}) OR (CMAKE_Fortran_COMPILER_ID STREQUAL GNU))
  set(_comp ${arith})
  if(NOT mpi)
    list(APPEND _comp mpiseq)
  endif()
  if(scotch)
    list(APPEND _comp Scotch)
  endif()
  if(metis)
    list(APPEND _comp METIS)
  endif()
  if(OpenMP_FOUND)
    list(APPEND _comp OpenMP)
  endif()

  if(autobuild)
    find_package(MUMPS COMPONENTS ${_comp})
  else()
    find_package(MUMPS COMPONENTS ${_comp} REQUIRED)
  endif()
endif()

if(NOT MUMPS_FOUND)
  set(mumps_external true CACHE BOOL "autobuild Mumps")

  # necessary since CMAKE_ARGS is broken in general
  set(parallel ${mpi} CACHE BOOL "Mumps parallel == Gemini mpi")
  set(MUMPS_BUILD_TESTING false CACHE BOOL "mumps disable tests")

  include(FetchContent)

  FetchContent_Declare(MUMPS_proj
    GIT_REPOSITORY https://github.com/scivision/mumps.git
    GIT_TAG v5.3.5.0
    CMAKE_ARGS -Darith=${arith} -Dmetis:BOOL=${metis} -Dscotch:BOOL=${scotch} -Dopenmp:BOOL=false
  )

  FetchContent_MakeAvailable(MUMPS_proj)
endif()

target_link_libraries(MUMPS::MUMPS INTERFACE SCALAPACK::SCALAPACK LAPACK::LAPACK)
