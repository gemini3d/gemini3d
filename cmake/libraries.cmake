# MKL Lapack / Blas
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor

get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)

if(USEMKL OR DEFINED ENV{MKLROOT} OR CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  if(NOT C IN_LIST languages)
    enable_language(C)
  endif()

  find_package(LAPACK REQUIRED COMPONENTS IntelPar)
else()
  find_package(LAPACK REQUIRED COMPONENTS Netlib)
endif()

# MPI
find_package(MPI REQUIRED COMPONENTS Fortran)

