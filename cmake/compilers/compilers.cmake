include(CheckSourceCompiles)
include(CheckCompilerFlag)

# check C and Fortran compiler ABI compatibility

if(NOT abi_ok)
  message(CHECK_START "checking that C and Fortran compilers can link")
  try_compile(abi_ok ${CMAKE_CURRENT_BINARY_DIR}/abi_check ${CMAKE_CURRENT_LIST_DIR}/abi_check abi_check)
  if(abi_ok)
    message(CHECK_PASS "OK")
  else()
    message(FATAL_ERROR "ABI-incompatible: C compiler ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION} and Fortran compiler ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
  endif()
endif()


set(CMAKE_EXPORT_COMPILE_COMMANDS on)

# === check that the compiler has adequate Fortran 2008 support
# this is to mitigate confusing syntax error messages for new users

# clean out prior libs to avoid false fails
set(CMAKE_REQUIRED_LIBRARIES)
set(CMAKE_REQUIRED_INCLUDES)
set(CMAKE_REQUIRED_FLAGS)

if(dev)
check_source_compiles(Fortran
"program imp
implicit none (type, external)
end program"
f2018impnone
)
if(NOT f2018impnone)
  message(FATAL_ERROR "Compiler does not support Fortran 2018 IMPLICIT NONE (type, external): ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
endif()

check_source_compiles(Fortran
"program es2018
character :: x
error stop x
end program"
f2018errorstop
)
if(NOT f2018errorstop)
  message(FATAL_ERROR "Compiler does not support Fortran 2018 error stop with character variable: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
endif()
endif()

if(dev)
check_source_compiles(Fortran
"program f18_assumed_rank
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
end program"
f2018assumed_rank
)
endif(dev)
