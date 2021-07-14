add_library(neutral OBJECT neutral.f90 atmos.f90 interp.f90 perturb.f90 proj.f90 wind.f90)
target_compile_options(neutral PRIVATE ${gcc10opts})
target_link_libraries(neutral PRIVATE config const reader grid meshobj interp mpimod
hwm_ifc msis_ifc
timeutils MPI::MPI_Fortran)

# self-tests
add_executable(test_proj test_proj.f90)
target_link_libraries(test_proj PRIVATE pathlib const newton meshobj meshobj_dipole spherical geomagnetic neutral config const reader grid interp mpimod hwm_ifc msis_ifc timeutils autogrid MPI::MPI_Fortran)


