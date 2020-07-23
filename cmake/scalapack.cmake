# Finds Scalapack, tests, and if not found or broken, autobuild scalapack
if(scalapack_external)
  include(${CMAKE_CURRENT_LIST_DIR}/scalapack_external.cmake)
  return()
endif()

if(autobuild)
  find_package(SCALAPACK)
else()
  find_package(SCALAPACK REQUIRED)
endif()

if(SCALAPACK_FOUND)
  set(scalapack_external false CACHE BOOL "autobuild Scalapack")
else()
  include(${CMAKE_CURRENT_LIST_DIR}/scalapack_external.cmake)
  return()
endif()

if(lapack_external)
# can't run prebuild test with external libraries not yet built.
  return()
endif()

# -- verify links

set(CMAKE_REQUIRED_FLAGS)
set(CMAKE_REQUIRED_LINK_OPTIONS)
set(CMAKE_REQUIRED_INCLUDES)
set(CMAKE_REQUIRED_LIBRARIES SCALAPACK::SCALAPACK LAPACK::LAPACK MPI::MPI_Fortran)
# MPI needed for ifort
include(CheckFortranSourceCompiles)

if("s" IN_LIST arith)

set(_code "
program test_real32
use, intrinsic :: iso_fortran_env, only: real32
implicit none
integer :: ictxt, myid, nprocs
real(real32) :: eps
real(real32), external :: pslamch
external :: blacs_pinfo, blacs_get, blacs_gridinit, blacs_gridexit, blacs_exit

call blacs_pinfo(myid, nprocs)
call blacs_get(-1, 0, ictxt)
call blacs_gridinit(ictxt, 'C', nprocs, 1)
eps = pslamch(ictxt, 'E')
call blacs_gridexit(ictxt)
call blacs_exit(0)

end program")

check_fortran_source_compiles(${_code} SCALAPACK_real32_link SRC_EXT f90)

if(NOT SCALAPACK_real32_link)
  message(STATUS "Scalapack ${SCALAPACK_LIBRARIES} not building with LAPACK ${LAPACK_LIBRARIES} and ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
  if(NOT autobuild)
    message(FATAL_ERROR "autobuild=off, so cannot proceed")
  endif()
  include(${CMAKE_CURRENT_LIST_DIR}/scalapack_external.cmake)
  return()
endif()
endif()

if("d" IN_LIST arith)

set(_code "
program test_real64
use, intrinsic :: iso_fortran_env, only: real64
implicit none
integer :: ictxt, myid, nprocs
real(real64) :: eps
real(real64), external :: pdlamch
external :: blacs_pinfo, blacs_get, blacs_gridinit, blacs_gridexit, blacs_exit

call blacs_pinfo(myid, nprocs)
call blacs_get(-1, 0, ictxt)
call blacs_gridinit(ictxt, 'C', nprocs, 1)
eps = pdlamch(ictxt, 'E')
call blacs_gridexit(ictxt)
call blacs_exit(0)

end program")

check_fortran_source_compiles(${_code} SCALAPACK_real64_link SRC_EXT f90)

if(NOT SCALAPACK_real64_link)
  message(STATUS "Scalapack ${SCALAPACK_LIBRARIES} not building with LAPACK ${LAPACK_LIBRARIES} and ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
  if(NOT autobuild)
    message(FATAL_ERROR "autobuild=off, so cannot proceed")
  endif()
  include(${CMAKE_CURRENT_LIST_DIR}/scalapack_external.cmake)
  return()
endif()
endif()


# FIXME: figure out SCALAPACK_run with non-default libgfortran.so
# as this test fails in that case. Is it an issue with LD_LIBRARY_PATH?
# if so, may be able to pass via CMAKE_REQUIRED_LINK_OPTIONS of check_fortran_source_runs()
# to underlying try_run() link_options. This would have to be tested on a system where
# this check currently fails, such as one using a non-default Gfortran.
include(CheckFortranSourceRuns)
check_fortran_source_runs(${_code} SCALAPACK_run)
if(NOT SCALAPACK_run)
  message(STATUS "Scalapack ${SCALAPACK_LIBRARIES} not running with LAPACK ${LAPACK_LIBRARIES} and ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
endif()
