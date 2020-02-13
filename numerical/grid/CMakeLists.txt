add_library(mesh OBJECT mesh.f90)
target_link_libraries(mesh PRIVATE const)

add_library(grid OBJECT grid.f90 read.f90)
target_link_libraries(grid PUBLIC const reader mesh mpimod MPI::MPI_Fortran) # PUBLIC necessary for ifort, even with MPI target
if(hdf5)
  target_sources(grid PRIVATE readgrid_hdf5.f90)
  target_link_libraries(grid PRIVATE h5fortran::h5fortran)
elseif(netcdf)
  target_sources(grid PRIVATE readgrid_nc4.f90)
  target_link_libraries(grid PRIVATE nc4fortran::nc4fortran)
else()
  target_sources(grid PRIVATE readgrid_raw.f90)
endif()
