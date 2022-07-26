add_library(grid OBJECT grid.f90 readgrid_hdf5.f90)
target_link_libraries(grid PRIVATE const reader meshobj
h5fortran::h5fortran
)

add_library(grid_mpi OBJECT grid_mpi.f90 read_mpi.f90 check_mpi.f90
$<TARGET_OBJECTS:grid>
)
target_compile_options(grid_mpi PRIVATE ${gfortran_opts})
target_link_libraries(grid_mpi PRIVATE const meshobj meshobj_cart meshobj_dipole mpimod
h5fortran::h5fortran
)
