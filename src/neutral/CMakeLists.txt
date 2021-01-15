add_library(neutral OBJECT neutral.f90 atmos.f90 interp.f90 perturb.f90 proj.f90)
target_compile_options(neutral PRIVATE ${gcc10opts})
target_link_libraries(neutral PRIVATE const reader grid mesh interp mpimod msis_ifc timeutils MPI::MPI_Fortran)
