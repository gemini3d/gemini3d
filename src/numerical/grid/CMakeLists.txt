add_library(grid OBJECT grid.f90 readgrid_hdf5.f90 readgrid_nc4.f90 readgrid_raw.f90
$<$<NOT:$<BOOL:${netcdf}>>:$<TARGET_OBJECTS:nc4fortran_dummy>>
)
target_link_libraries(grid PRIVATE const reader meshobj 
h5fortran::h5fortran
$<$<BOOL:${netcdf}>:nc4fortran::nc4fortran>
$<$<BOOL:${HDF5_HAVE_PARALLEL}>:MPI::MPI_Fortran>
)

add_library(grid_mpi OBJECT grid_mpi.f90 read_mpi.f90 check_mpi.f90)
target_compile_options(grid_mpi PRIVATE ${gfortran_opts})
target_link_libraries(grid_mpi PRIVATE const grid meshobj meshobj_cart meshobj_dipole mpimod
h5fortran::h5fortran
$<$<BOOL:${netcdf}>:nc4fortran::nc4fortran>
$<$<BOOL:${HDF5_HAVE_PARALLEL}>:MPI::MPI_Fortran>
)
