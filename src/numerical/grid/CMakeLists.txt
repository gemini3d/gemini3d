add_library(grid OBJECT grid.f90 readgrid_hdf5.f90
grid_mpi.f90 read_mpi.f90 check_mpi.f90
)
target_link_libraries(grid PRIVATE const reader meshobj meshobj_cart meshobj_dipole
mpimod
h5fortran::h5fortran
MPI::MPI_Fortran
)
