#if(NOT CMAKE_Fortran_COMPILER_VERSION VERSION_EQUAL 7.2.0)
list(APPEND FFLAGS -march=native)
#endif()

list(APPEND FFLAGS -fimplicit-none)
list(APPEND FFLAGS -Wall -Wpedantic -Wextra)

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  list(APPEND FFLAGS -Werror=array-bounds -fcheck=all)
  # list(APPEND FFLAGS -ffpe-trap=invalid,zero,overflow)#,underflow)
else()
  list(APPEND FFLAGS -Wno-unused-dummy-argument -Wno-unused-variable -Wno-unused-function)
endif()

# enforce Fortran 2018 standard
if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 8)
   list(APPEND FFLAGS -std=f2018)
endif()