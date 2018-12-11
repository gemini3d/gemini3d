cmake_policy(SET CMP0074 NEW)

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  add_compile_options(-g -O0)
else()
  add_compile_options(-g -O3)
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  add_compile_options(-qopenmp)  # undefined reference to `omp_get_max_threads'
  add_link_options(-parallel) # undefined reference to `__kmpc_begin'


  if(CMAKE_BUILD_TYPE STREQUAL Debug)
    add_compile_options(-debug extended -check all -heap-arrays -fp-stack-check)
  endif()
  add_compile_options(-warn nounused -traceback -stand f08 
    -diag-disable 5268)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  # -fdefault-real-8  after literals are fixed to "e" or "wp"
  add_compile_options(-march=native -fimplicit-none)
  add_compile_options(-Wall -Wpedantic -Wextra)
  
  if(CMAKE_BUILD_TYPE STREQUAL Debug)
    add_compile_options(-fcheck=all)
    # add_compile_options(-ffpe-trap=invalid,zero,overflow)#,underflow)
  else()
    add_compile_options(-Wno-unused-dummy-argument -Wno-unused-variable -Wno-unused-function)
  endif()
  
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 6)
     add_compile_options(-std=f2008)
  elseif(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 8)
     add_compile_options(-std=f2018)
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL PGI)

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Cray)

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL XL)

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Flang) 
  add_compile_options(-Mallocatable=03)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL NAG)
  add_compile_options(-u -C=all -f2008)
endif()
