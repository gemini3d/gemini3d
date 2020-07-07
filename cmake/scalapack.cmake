# Finds Scalapack, tests, and if not found or broken, autobuild scalapack
if(scalapack_external)
  include(${CMAKE_CURRENT_LIST_DIR}/scalapack_external.cmake)
  return()
endif()

find_package(SCALAPACK)

if(NOT SCALAPACK_FOUND)
  include(${CMAKE_CURRENT_LIST_DIR}/scalapack_external.cmake)
  set(scalapack_external true CACHE BOOL "autobuild Scalapack")
  return()
else()
  set(scalapack_external false CACHE BOOL "autobuild Scalapack")
endif()

if(lapack_external)
# can't run prebuild test with external libraries not yet built.
  return()
endif()

# -- verify Scalapack links

set(CMAKE_REQUIRED_LIBRARIES SCALAPACK::SCALAPACK LAPACK::LAPACK MPI::MPI_Fortran)
# MPI needed for ifort

if("s" IN_LIST arith)
  set(_code "use, intrinsic :: iso_fortran_env, only: real32
implicit none
integer :: ictxt, myid, nprocs
real(real32) :: eps
real(real32), external :: pslamch

call blacs_pinfo(myid, nprocs)
call blacs_get(-1, 0, ictxt)
call blacs_gridinit(ictxt, 'C', nprocs, 1)
eps = pslamch(ictxt, 'E')
call blacs_gridexit(ictxt)
call blacs_exit(0)

end program")
elseif("d" IN_LIST arith)
  set(_code "use, intrinsic :: iso_fortran_env, only: real64
implicit none
integer :: ictxt, myid, nprocs
real(real64) :: eps
real(real64), external :: pdlamch

call blacs_pinfo(myid, nprocs)
call blacs_get(-1, 0, ictxt)
call blacs_gridinit(ictxt, 'C', nprocs, 1)
eps = pdlamch(ictxt, 'E')
call blacs_gridexit(ictxt)
call blacs_exit(0)

end program")
else()
  message(STATUS "SKIP: Scalapack test arith: ${arith}")
  return()
endif()
include(CheckFortranSourceCompiles)
check_fortran_source_compiles(${_code} SCALAPACK_Compiles_OK SRC_EXT f90)
if(NOT SCALAPACK_Compiles_OK)
  message(STATUS "Scalapack ${SCALAPACK_LIBRARIES} not building with LAPACK ${LAPACK_LIBRARIES} and ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
  include(${CMAKE_CURRENT_LIST_DIR}/scalapack_external.cmake)
  set(scalapack_external true CACHE BOOL "autobuild Scalapack")
endif()

# FIXME: figure out SCALAPACK_Runs_OK with non-default libgfortran.so
# as this test fails in that case. Is it an issue with LD_LIBRARY_PATH?
# if so, may be able to pass via CMAKE_REQUIRED_LINK_OPTIONS of check_fortran_source_runs()
# to underlying try_run() link_options. This would have to be tested on a system where
# this check currently fails, such as one using a non-default Gfortran.
include(CheckFortranSourceRuns)
check_fortran_source_runs(${_code} SCALAPACK_Runs_OK SRC_EXT f90)
if(NOT SCALAPACK_Runs_OK)
  message(STATUS "Scalapack ${SCALAPACK_LIBRARIES} not running with LAPACK ${LAPACK_LIBRARIES} and ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
endif()
