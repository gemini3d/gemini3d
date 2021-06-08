add_library(grid grid.f90 read.f90 check.f90
readgrid_hdf5.f90 readgrid_nc4.f90 readgrid_raw.f90)
# not OBJECT to avoid problems with CMake 3.13
target_compile_options(grid PRIVATE ${gcc10opts})
target_link_libraries(grid PRIVATE const reader meshobj meshobj_cart meshobj_dipole mpimod MPI::MPI_Fortran
h5fortran::h5fortran nc4fortran::nc4fortran)
