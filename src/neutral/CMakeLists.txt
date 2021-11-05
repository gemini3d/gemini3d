add_library(neutral OBJECT neutral.f90 atmos.f90 interp.f90 perturb.f90 proj.f90 wind.f90)
target_compile_options(neutral PRIVATE ${gfortran_opts})
target_link_libraries(neutral PRIVATE config const reader pathlib grid meshobj interp mpimod
hwm_ifc msis_ifc
timeutils
neutraldata neutraldata3D
h5fortran::h5fortran
MPI::MPI_Fortran
)

# self-tests
add_executable(test_proj test_proj.f90
$<TARGET_OBJECTS:newton> $<TARGET_OBJECTS:spherical> $<TARGET_OBJECTS:geomagnetic>
$<TARGET_OBJECTS:meshobj> $<TARGET_OBJECTS:meshobj_cart> $<TARGET_OBJECTS:meshobj_dipole>
$<TARGET_OBJECTS:grid> $<TARGET_OBJECTS:autogrid>
$<TARGET_OBJECTS:neutral> $<TARGET_OBJECTS:interp>
)
target_include_directories(test_proj PRIVATE ${PROJECT_BINARY_DIR}/include)
target_link_libraries(test_proj PRIVATE config const reader
timeutils pathlib mpimod const
hwm_ifc
msis_ifc
MPI::MPI_Fortran
)
