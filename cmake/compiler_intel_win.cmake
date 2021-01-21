add_compile_options(/QxHost)  # like -march=native

# NOTE: /nologo /libs:dll are set by CMake in cmake_fortran_flags
string(APPEND CMAKE_Fortran_FLAGS " /warn:declarations /traceback")
string(APPEND CMAKE_Fortran_FLAGS " /Qopenmp")
string(APPEND CMAKE_Fortran_FLAGS " /heap-arrays")  # necessary for stack overflow avoid, even for 2D
string(APPEND CMAKE_Fortran_FLAGS_DEBUG " /check:bounds")
string(APPEND CMAKE_Fortran_FLAGS " /warn:nounused /Qdiag-disable:5268 /Qdiag-disable:7712")
#add_link_options(/Qparallel)

# enforce Fortran 2018 standard
# string(APPEND CMAKE_Fortran_FLAGS " /stand:f18")  # too many false warnings
