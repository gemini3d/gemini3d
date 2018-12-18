# MKL Lapack / Blas
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor

get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)

if(DEFINED ENV{MKLROOT})
  if(NOT C IN_LIST languages)
    message(FATAL_ERROR "MKL/Intel requires project to also have C language enabled")
  endif()
#  add_link_options(-Wl,--no-as-needed)

  set(BLA_VENDOR Intel10_64lp_seq)
  #include_directories($ENV{MKLROOT}/include)
endif()


# Lapack95 - do this BEFORE regular Lapack!
set(BLA_F95 ON)
find_package(LAPACK QUIET)  # sets LAPACK95_FOUND, LAPACK95_LIBRARIES

# Lapack
set(BLA_F95 OFF)
#find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# MPI
find_package(MPI REQUIRED COMPONENTS Fortran)

