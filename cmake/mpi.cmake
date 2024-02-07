include(CheckSourceCompiles)

if(NOT DEFINED MPI_ROOT AND DEFINED ENV{MPI_ROOT})
  set(MPI_ROOT $ENV{MPI_ROOT})
endif()
if(MPI_ROOT)
  message(STATUS "Using MPI_ROOT=${MPI_ROOT}")
endif()

set(MPI_DETERMINE_LIBRARY_VERSION true)

find_package(MPI COMPONENTS C CXX Fortran REQUIRED)

message(STATUS "${MPI_Fortran_LIBRARY_VERSION_STRING}")
message(STATUS "MPI libs: ${MPI_Fortran_LIBRARIES}")
message(STATUS "MPI include: ${MPI_Fortran_INCLUDE_DIRS}")
message(STATUS "MPI compile flags: ${MPI_Fortran_COMPILER_FLAGS}")
message(STATUS "MPI link flags: ${MPI_Fortran_LINK_FLAGS}")

include(${CMAKE_CURRENT_LIST_DIR}/openmpi.cmake)

if(MPI_Fortran_HAVE_F08_MODULE)
  return()
endif()

set(CMAKE_REQUIRED_LIBRARIES MPI::MPI_Fortran)

# sometimes factory FindMPI.cmake doesn't define this
message(CHECK_START "Checking for Fortran MPI-3 binding")
check_source_compiles(Fortran
"program test
use mpi_f08, only : mpi_comm_rank, mpi_real, mpi_comm_world, mpi_init, mpi_finalize
implicit none
call mpi_init
call mpi_finalize
end program"
MPI_Fortran_HAVE_F08_MODULE
)

if(MPI_Fortran_HAVE_F08_MODULE)
  message(CHECK_PASS "yes")
else()
  message(CHECK_FAIL "no -- fallback to MPI-2 interface")
endif()
