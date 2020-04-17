include(CheckFortranSourceCompiles)
include(CheckFortranSourceRuns)
include(CheckFortranCompilerFlag)

# === check that the compiler has adequate Fortran 2008 support
# this is to mitigate confusing syntax error messages for new users

check_fortran_source_compiles("implicit none (external); end" f2018impnone SRC_EXT f90)
if(NOT f2018impnone)
  message(FATAL_ERROR "Compiler does not support Fortran 2018 IMPLICIT NONE (EXTERNAL): ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
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

# Do these before compiler options so options don't goof up finding
# === OpenMP
# optional, for possible MUMPS speedup
if(openmp)
  find_package(OpenMP COMPONENTS C Fortran)
endif()

# === MPI
# MPI is used throughout Gemini
find_package(MPI REQUIRED COMPONENTS Fortran)
# sometimes a mismatch exists between CMake's found MPI and the compiler.
# Catch this here to avoid confusing build errors.
# Runtime errors will be caught by test_mpi in ctest.
set(CMAKE_REQUIRED_INCLUDES ${MPI_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES MPI::MPI_Fortran)
check_fortran_source_compiles("use mpi; end" MPI_OK SRC_EXT f90)
if(NOT MPI_OK)
  message(FATAL_ERROR "${MPI_Fortran_LIBRARIES} not working with ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
endif()
