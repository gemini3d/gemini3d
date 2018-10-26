## ----- release or debug?

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  add_compile_options(-g -O0)
else()
  add_compile_options(-g -O3)
endif()

## ------- Compiler options
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
  # -r8  after literals are fixed to "e" or "wp"
  if(CMAKE_BUILD_TYPE STREQUAL Debug)
     add_compile_options(-debug extended -check all -heap-arrays)
  endif()
  add_compile_options(-warn nounused -traceback -fp-stack-check)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
  # -fdefault-real-8  after literals are fixed to "e" or "wp"
  add_compile_options(-mtune=native -fimplicit-none)
  
  if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 6)
    add_compile_options(-Wall -Wpedantic -Wextra -fexceptions -fbacktrace -Wno-unused-dummy-argument -Wno-unused-variable -Wno-unused-function)
  endif()
# -fstack-protector
# -ffpe-trap=invalid,zero,overflow)#,underflow)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL PGI)
  add_compile_options(-Mallocatable=03)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Flang) 
  add_compile_options(-Mallocatable=03)
  link_libraries(-static-flang-libs)
endif()
