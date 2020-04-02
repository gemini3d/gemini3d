include(CheckFortranSourceCompiles)
include(CheckFortranSourceRuns)

# === check that the compiler has adequate Fortran 2008 support
# this is to mitigate confusing syntax error messages for new users

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

# === compiler setup
# feel free to add more compiler_*.cmake
if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  include(${CMAKE_CURRENT_LIST_DIR}/compiler_intel.cmake)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  include(${CMAKE_CURRENT_LIST_DIR}/compiler_gnu.cmake)
endif()
