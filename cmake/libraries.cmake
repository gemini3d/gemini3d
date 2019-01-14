# MKL Lapack / Blas
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor

get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)

if(DEFINED ENV{MKLROOT})
  if(NOT C IN_LIST languages)
    message(FATAL_ERROR "MKL/Intel requires project to also have C language enabled")
  endif()
#  add_link_options(-Wl,--no-as-needed)

  #include_directories($ENV{MKLROOT}/include)
  
  find_package(LAPACK REQUIRED COMPONENTS IntelPar)
else()
  find_package(LAPACK REQUIRED)
endif()

# MPI
find_package(MPI REQUIRED COMPONENTS Fortran)

