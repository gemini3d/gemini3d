include(CheckFortranCompilerFlag)

# --- abi check: C++ and Fortran compiler ABI compatibility

function(abi_check)
if(NOT abi_compile)
  message(CHECK_START "checking that C, C++, and Fortran compilers can link")
  try_compile(abi_compile
  ${CMAKE_CURRENT_BINARY_DIR}/abi_check ${CMAKE_CURRENT_LIST_DIR}/abi_check
  abi_check
  OUTPUT_VARIABLE abi_output
  )
  if(abi_output MATCHES "ld: warning: could not create compact unwind for")
    message(WARNING "C++ exception handling will not work reliably due to incompatible compilers:
    C++ compiler ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}
    Fortran compiler ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}"
    )
  endif()

  if(abi_compile)
    message(CHECK_PASS "OK")
  else()
    message(FATAL_ERROR "ABI-incompatible compilers:
    C compiler ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}
    C++ compiler ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}
    Fortran compiler ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}"
    )
  endif()
endif()
endfunction(abi_check)
abi_check()

if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  include(${CMAKE_CURRENT_LIST_DIR}/intel.cmake)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  include(${CMAKE_CURRENT_LIST_DIR}/gnu.cmake)
endif()
