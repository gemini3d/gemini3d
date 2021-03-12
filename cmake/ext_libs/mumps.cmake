# MUMPS -- use cmake -DMUMPS_ROOT= for hint.
#
# Intel MKL-compiled MUMPS requires at the linker for the main executable:
# mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 mkl_intel_lp64 mkl_intel_thread mkl_core
#
# easily obtain MUMPS without compiling:
# CentOS 6/7 EPEL: yum install mumps-devel
# Ubuntu / Debian: apt install libmumps-dev

include(FetchContent)

# --- prereqs
include(${CMAKE_CURRENT_LIST_DIR}/lapack.cmake)

if(mpi)
  include(${CMAKE_CURRENT_LIST_DIR}/scalapack.cmake)

  find_package(HWLOC)  # need here for cmake/summary.cmake scope
endif()
# --- MUMPS

if(NOT mumps_external AND (MUMPS_ROOT OR (DEFINED ENV{MUMPS_ROOT}) OR (CMAKE_Fortran_COMPILER_ID STREQUAL GNU)))
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

if(MUMPS_FOUND OR TARGET MUMPS::MUMPS)
  return()
endif()

set(mumps_external true CACHE BOOL "build Mumps")

# cache to override mumps option(parallel on)
set(parallel ${mpi} CACHE BOOL "Mumps parallel = Gemini mpi")

FetchContent_Declare(MUMPS
  GIT_REPOSITORY ${mumps_git}
  GIT_TAG ${mumps_tag}
  CMAKE_ARGS -Darith=${arith} -Dmetis:BOOL=${metis} -Dscotch:BOOL=${scotch} -Dopenmp:BOOL=false)

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.14)
  FetchContent_MakeAvailable(MUMPS)
elseif(NOT mumps_POPULATED)
  FetchContent_Populate(MUMPS)
  add_subdirectory(${mumps_SOURCE_DIR} ${mumps_BINARY_DIR})
endif()
