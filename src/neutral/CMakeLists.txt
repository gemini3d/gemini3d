add_library(neutral OBJECT neutral.f90 atmos.f90)
target_link_libraries(neutral PRIVATE const reader grid mesh interp mpimod msis timeutils MPI::MPI_Fortran)
