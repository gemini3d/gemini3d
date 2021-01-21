add_compile_options(-xHost)  # like -march=native

string(APPEND CMAKE_Fortran_FLAGS " -warn declarations -traceback")  # -warn all or -warn gets mixed with -qopenmp with CMake 3.14.2
string(APPEND CMAKE_Fortran_FLAGS " -qopenmp")  # undefined reference to `omp_get_max_threads'
string(APPEND CMAKE_Fortran_FLAGS " -heap-arrays")  # (is this needed on Linux?) stack overflow avoid
#string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -check all")
#string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -debug extended -check all -fpe0 -fp-stack-check")
string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -check bounds")
string(APPEND CMAKE_Fortran_FLAGS " -warn nounused -diag-disable 5268 -diag-disable 7712")
add_link_options(-parallel) # without this, error: undefined reference to `__kmpc_begin'

# enforce Fortran 2018 standard
string(APPEND CMAKE_Fortran_FLAGS " -stand f18")
