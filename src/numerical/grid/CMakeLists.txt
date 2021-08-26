add_library(grid OBJECT grid.f90 read.f90 check.f90
readgrid_hdf5.f90 readgrid_nc4.f90 readgrid_raw.f90)
target_compile_options(grid PRIVATE ${gfortran_opts})
target_link_libraries(grid PRIVATE const reader meshobj meshobj_cart meshobj_dipole mpimod
h5fortran::h5fortran nc4fortran::nc4fortran
MPI::MPI_Fortran)
