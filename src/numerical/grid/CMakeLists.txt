configure_file(mesh.in.f90 mesh.f90)
add_library(mesh OBJECT ${CMAKE_CURRENT_BINARY_DIR}/mesh.f90)

add_library(grid OBJECT grid.f90 read.f90 check.f90
readgrid_hdf5.f90 readgrid_nc4.f90 readgrid_raw.f90)
target_compile_options(grid PRIVATE ${gcc10opts})
 # PUBLIC necessary for ifort, even with MPI target
target_link_libraries(grid PRIVATE const reader mesh mpimod MPI::MPI_Fortran
h5fortran::h5fortran nc4fortran::nc4fortran)
