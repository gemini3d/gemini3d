# MKL Lapack / Blas
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
# if not finding Intel/MKL, be sure C language is enabled in CMakeLists.txt
if(USEMKL OR CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  if(NOT DEFINED ENV{MKLROOT})
    message(FATAL_ERROR "MKLROOT must be defined")
  endif()
  
  add_link_options(-Wl,--no-as-needed)

  set(BLA_VENDOR Intel10_64lp_seq)
  include_directories($ENV{MKLROOT}/include)
endif()


# Lapack95 - do this BEFORE regular Lapack!
if(LIB_DIR AND NOT USEMKL AND NOT CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  list(APPEND LAPACK95_ROOT ${LIB_DIR}/LAPACK95)
endif()
set(BLA_F95 ON)
find_package(LAPACK QUIET)  # sets LAPACK95_FOUND, LAPACK95_LIBRARIES

# Lapack
set(BLA_F95 OFF)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)


# MPI
find_package(MPI REQUIRED COMPONENTS Fortran)

# SCALAPACK
if(LIB_DIR AND NOT USEMKL AND NOT CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  list(APPEND SCALAPACK_ROOT ${LIB_DIR}/scalapack)
endif()
find_package(SCALAPACK REQUIRED)


