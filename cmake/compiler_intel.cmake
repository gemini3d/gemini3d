
if(WIN32)
  # workaround: /nologo /libs:dll are set by CMake in cmake_fortran_flags
  string(APPEND CMAKE_Fortran_FLAGS " /warn:declarations /traceback")
  string(APPEND CMAKE_Fortran_FLAGS " /Qopenmp")
  string(APPEND CMAKE_Fortran_FLAGS " /heap-arrays")  # necessary for stack overflow avoid, even for 2D
  string(APPEND CMAKE_Fortran_FLAGS_DEBUG " /check:bounds")
  string(APPEND CMAKE_Fortran_FLAGS " /warn:nounused /Qdiag-disable:5268 /Qdiag-disable:7712")
  #add_link_options(/Qparallel)
else()
  string(APPEND CMAKE_Fortran_FLAGS " -warn declarations -traceback")  # -warn all or -warn gets mixed with -qopenmp with CMake 3.14.2
  string(APPEND CMAKE_Fortran_FLAGS " -qopenmp")  # undefined reference to `omp_get_max_threads'
  string(APPEND CMAKE_Fortran_FLAGS " -heap-arrays")  # (is this needed on Linux?) stack overflow avoid
  #string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -check all")
  #string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -debug extended -check all -fpe0 -fp-stack-check")
  string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -check bounds")
  string(APPEND CMAKE_Fortran_FLAGS " -warn nounused -diag-disable 5268 -diag-disable 7712")
  add_link_options(-parallel) # undefined reference to `__kmpc_begin'
endif(WIN32)


# enforce Fortran 2018 standard
if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 19)
  if(WIN32)
    string(APPEND CMAKE_Fortran_FLAGS " /stand:f18")
  else()
    string(APPEND CMAKE_Fortran_FLAGS " -stand f18")
  endif()
endif()
