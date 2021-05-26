include(CheckFortranSourceCompiles)
include(CheckFortranCompilerFlag)

# check C and Fortran compiler ABI compatibility

if(NOT abi_ok)
  message(CHECK_START "checking that C and Fortran compilers can link")
  try_compile(abi_ok ${CMAKE_CURRENT_BINARY_DIR}/abi_check ${CMAKE_CURRENT_LIST_DIR}/abi_check abi_check)
  if(abi_ok)
    message(CHECK_PASS "OK")
  else()
    message(FATAL ERROR "C compiler {CMAKE_C_COMPILER_ID} {CMAKE_C_COMPILER_VERSION} and Fortran compiler ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION} are ABI-incompatible.")
  endif()
endif()



set(CMAKE_EXPORT_COMPILE_COMMANDS on)
set(CMAKE_CONFIGURATION_TYPES "Release;Debug" CACHE STRING "Build type selections")


# === check that the compiler has adequate Fortran 2008 support
# this is to mitigate confusing syntax error messages for new users

# clean out prior libs to avoid false fails
set(CMAKE_REQUIRED_LIBRARIES)
set(CMAKE_REQUIRED_INCLUDES)
set(CMAKE_REQUIRED_FLAGS)

check_fortran_source_compiles("implicit none (type, external); end" f2018impnone SRC_EXT f90)
if(NOT f2018impnone)
  message(FATAL_ERROR "Compiler does not support Fortran 2018 IMPLICIT NONE (type, external): ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
endif()

check_fortran_source_compiles("program es2018
character :: x
error stop x
end program" f2018errorstop SRC_EXT f90)
if(NOT f2018errorstop)
  message(FATAL_ERROR "Compiler does not support Fortran 2018 error stop with character variable: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
endif()

check_fortran_source_compiles("program f18_assumed_rank
implicit none (type, external)
contains
subroutine ranker(A)
integer, intent(in) :: A(..)
select rank(A)
  rank (0)
    print *, rank(A)
  rank default
    print *, rank(A)
end select
end subroutine ranker
end program" f2018assumed_rank SRC_EXT f90)

# --- MSISE00 and MSIS 2.0 require legacy workaround due to non-standard Fortran code

if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  # Gfortran >= 8 need -Wno-pedantic to allow mismatched array size inhernet to MSIS.
  # "-w" doesn't disable pedantic
  set(msis_flags -w -std=legacy -Wno-pedantic -fno-implicit-none -Wno-error=array-bounds -fcheck=no-all)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Intel OR CMAKE_Fortran_COMPILER_ID STREQUAL IntelLLVM)
  set(msis_flags -nowarn)
endif()

# Do these before compiler options so options don't goof up finding

if(mpi)
  include(${PROJECT_SOURCE_DIR}/cmake/ext_libs/mpi.cmake)
endif(mpi)

if(NOT mpi)
  add_subdirectory(src/vendor/mpi_stubs)
endif(NOT mpi)

# optional, for possible MUMPS speedup
if(openmp)
find_package(OpenMP COMPONENTS C Fortran)
endif()
