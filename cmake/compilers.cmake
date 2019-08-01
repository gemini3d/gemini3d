cmake_policy(SET CMP0074 NEW)
cmake_policy(SET CMP0075 NEW)
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
  cmake_policy(SET CMP0076 NEW)
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)

  if(WIN32)
    set(FFLAGS /warn:declarations /traceback)
    list(APPEND FFLAGS /Qopenmp)
  else()
    set(FFLAGS -warn declarations -traceback)  # -warn all or -warn gets mixed with -qopenmp with CMake 3.14.2
    list(APPEND FFLAGS -qopenmp)  # undefined reference to `omp_get_max_threads'
  endif()

  cmake_policy(VERSION 3.13)
  if(WIN32)
    #add_link_options(/Qparallel)
  else()
    add_link_options(-parallel) # undefined reference to `__kmpc_begin'
  endif()

  if(CMAKE_BUILD_TYPE STREQUAL Debug)
    if(WIN32)
      list(APPEND FFLAGS /check:bounds)
      list(APPEND FFLAGS /heap-arrays)  # stack overflow even on tiny "Test2D" case without this option
    else()
      #list(APPEND FFLAGS -check all)
      #list(APPEND FFLAGS -debug extended -check all -heap-arrays -fpe0 -fp-stack-check)
      list(APPEND FFLAGS -check bounds)
    endif()
  endif()

  if(WIN32)
    list(APPEND FFLAGS /warn:nounused /Qdiag-disable:5268)
  else()
    list(APPEND FFLAGS -warn nounused -diag-disable 5268)
  endif()

  if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 19)
    if(WIN32)
      list(APPEND FFLAGS /stand:f18)
    else()
      list(APPEND FFLAGS -stand f18)
    endif()
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  # -fdefault-real-8  after literals are fixed to "e" or "wp"
  if(NOT CMAKE_Fortran_COMPILER_VERSION VERSION_EQUAL 7.2.0)
    list(APPEND FFLAGS -march=native)
  endif()

  list(APPEND FFLAGS -fimplicit-none)
  list(APPEND FFLAGS -Wall -Wpedantic -Wextra)

  if(CMAKE_BUILD_TYPE STREQUAL Debug)
    list(APPEND FFLAGS -fcheck=all)
    # list(APPEND FFLAGS -ffpe-trap=invalid,zero,overflow)#,underflow)
  else()
    list(APPEND FFLAGS -Wno-unused-dummy-argument -Wno-unused-variable -Wno-unused-function)
  endif()

  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 8)
     list(APPEND FFLAGS -std=f2018)
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL PGI)

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Cray)

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL XL)

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Flang)

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL NAG)
  list(APPEND FFLAGS -u -C=all -f2008)
endif()
