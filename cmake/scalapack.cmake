# Finds Scalapack, tests, and if not found or broken, autobuild scalapack

include(${CMAKE_CURRENT_LIST_DIR}/lapack.cmake)

find_package(SCALAPACK)

set(scalapack_external false)
if(NOT SCALAPACK_FOUND)
  include(${CMAKE_CURRENT_LIST_DIR}/scalapack_external.cmake)
  set(scalapack_external true)
endif()

if(scalapack_external OR lapack_external)
# can't run prebuild test with external libraries not yet built.
  return()
endif()
# -- verify Scalapack links

set(CMAKE_REQUIRED_INCLUDES ${SCALAPACK_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${SCALAPACK_LIBRARIES} ${LAPACK_LIBRARIES} MPI::MPI_Fortran)
# MPI needed for ifort

# FIXME: make real32 test
file(READ ${CMAKE_SOURCE_DIR}/tests/test_scalapack_d.f90 _code)

check_fortran_source_compiles(${_code} SCALAPACK_Compiles_OK SRC_EXT f90)
if(NOT SCALAPACK_Compiles_OK)
  message(WARNING "Scalapack ${SCALAPACK_LIBRARIES} not building with LAPACK ${LAPACK_LIBRARIES} and ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
  include(${CMAKE_CURRENT_LIST_DIR}/scalapack_external.cmake)
  set(scalapack_external true)
endif()

# FIXME: figure out SCALAPACK_Runs_OK with non-default libgfortran.so
# as this test fails in that case. Is it an issue with LD_LIBRARY_PATH?
# if so, may be able to pass via CMAKE_REQUIRED_LINK_OPTIONS of check_fortran_source_runs()
# to underlying try_run() link_options. This would have to be tested on a system where
# this check currently fails, such as one using a non-default Gfortran.
check_fortran_source_runs(${_code} SCALAPACK_Runs_OK SRC_EXT f90)
if(NOT SCALAPACK_Runs_OK)
  message(STATUS "Scalapack ${SCALAPACK_LIBRARIES} not running with LAPACK ${LAPACK_LIBRARIES} and ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
endif()
