set(CMAKE_Fortran_FLAGS "-march=native -fimplicit-none -Wall -Wextra ")

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  string(APPEND CMAKE_Fortran_FLAGS "-Werror=array-bounds -fcheck=all ")
  # string(APPEND CMAKE_Fortran_FLAGS "-ffpe-trap=invalid,zero,overflow ")#,underflow)
else()
  string(APPEND CMAKE_Fortran_FLAGS "-Wno-unused-dummy-argument -Wno-unused-variable -Wno-unused-function ")
endif()

# enforce Fortran 2018 standard
if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 8)
   string(APPEND CMAKE_Fortran_FLAGS "-std=f2018 -Wpedantic ")
endif()