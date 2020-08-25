include(CheckFortranSourceCompiles)
include(CheckFortranSourceRuns)
include(CheckFortranCompilerFlag)

# === check that the compiler has adequate Fortran 2008 support
# this is to mitigate confusing syntax error messages for new users

# clean out prior libs to avoid false fails
set(CMAKE_REQUIRED_LIBRARIES)
set(CMAKE_REQURIED_INCLUDES)
set(CMAKE_REQUIRED_FLAGS)
check_fortran_source_compiles("implicit none (type, external); end" f2018impnone SRC_EXT f90)
if(NOT f2018impnone)
  message(FATAL_ERROR "Compiler does not support Fortran 2018 IMPLICIT NONE (type, external): ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
endif()

check_fortran_source_compiles("block; end block; end" f2008block SRC_EXT f90)
if(NOT f2008block)
  message(FATAL_ERROR "Compiler does not support Fortran 2008 BLOCK syntax: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
endif()

check_fortran_source_compiles("
module A
interface
module subroutine C
end subroutine
end interface
end module
submodule (A) B
contains
module procedure C
end procedure
end submodule
program D
end program" f2008submodule SRC_EXT f90)
if(NOT f2008submodule)
  message(FATAL_ERROR "Compiler does not support Fortran 2008 SUBMODULE syntax: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
endif()

if(NOT mpi)
  set(MPI_OK false)
  return()
endif()

# Do these before compiler options so options don't goof up finding
# === OpenMP
# optional, for possible MUMPS speedup
if(openmp)
  find_package(OpenMP COMPONENTS C Fortran)
endif()

# === MPI
include(${CMAKE_CURRENT_LIST_DIR}/mpi.cmake)
