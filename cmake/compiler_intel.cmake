if(WIN32)
  set(FFLAGS /warn:declarations /traceback)
  list(APPEND FFLAGS /Qopenmp)
  list(APPEND FFLAGS /heap-arrays)  # necessary for stack overflow avoid
else()
  set(FFLAGS -warn declarations -traceback)  # -warn all or -warn gets mixed with -qopenmp with CMake 3.14.2
  list(APPEND FFLAGS -qopenmp)  # undefined reference to `omp_get_max_threads'
  list(APPEND FFLAGS -heap-arrays)  # (is this needed on Linux?) stack overflow avoid
endif(WIN32)

if(WIN32)
  #add_link_options(/Qparallel)
else()
  add_link_options(-parallel) # undefined reference to `__kmpc_begin'
endif(WIN32)

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  if(WIN32)
    list(APPEND FFLAGS /check:bounds)
  else()
    #list(APPEND FFLAGS -check all)
    #list(APPEND FFLAGS -debug extended -check all -fpe0 -fp-stack-check)
    list(APPEND FFLAGS -check bounds)
  endif()
endif()

# reduce build-time warning verbosity
if(WIN32)
  list(APPEND FFLAGS /warn:nounused /Qdiag-disable:5268 /Qdiag-disable:7712)
else()
  list(APPEND FFLAGS -warn nounused -diag-disable 5268 -diag-disable 7712)
endif(WIN32)

# enforce Fortran 2018 standard
if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 19)
  if(WIN32)
    list(APPEND FFLAGS /stand:f18)
  else()
    list(APPEND FFLAGS -stand f18)
  endif()
endif()