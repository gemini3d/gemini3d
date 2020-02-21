add_compile_options(-march=native -Wall -Wextra)

string(APPEND CMAKE_Fortran_FLAGS "  -fimplicit-none")

# -Wpedantic makes too many false positives, through Gfortran 9

string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -Werror=array-bounds -fcheck=all")
  # string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -ffpe-trap=invalid,zero,overflow")#,underflow)
string(APPEND CMAKE_Fortran_FLAGS " -Wno-unused-dummy-argument -Wno-unused-variable -Wno-unused-function")

# enforce Fortran 2018 standard
if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 8)
  string(APPEND CMAKE_Fortran_FLAGS " -std=f2018")
endif()