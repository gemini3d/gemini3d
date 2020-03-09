# NOTE: don't use -march=native as GCC doesn't support all CPU arches with that option.
add_compile_options(-mtune=native)

# keep Wall and Wextra in cmake_fortran_flags so FetchContent packages can override
# and thereby avoid useless megabytes of other project warnings
string(APPEND CMAKE_Fortran_FLAGS " -Wall -Wextra -fimplicit-none")

# -Wpedantic makes too many false positives, through Gfortran 9

string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -Werror=array-bounds -fcheck=all")
  # string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -ffpe-trap=invalid,zero,overflow")#,underflow)
string(APPEND CMAKE_Fortran_FLAGS " -Wno-unused-dummy-argument -Wno-unused-variable -Wno-unused-function")

# enforce Fortran 2018 standard
if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 8 AND
   NOT CMAKE_Fortran_COMPILER_VERSION VERSION_EQUAL 9.2.0)
  string(APPEND CMAKE_Fortran_FLAGS " -std=f2018")
endif()
