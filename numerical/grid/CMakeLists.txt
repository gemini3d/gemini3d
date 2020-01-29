add_library(mesh mesh.f90)
target_link_libraries(mesh PRIVATE const)
target_include_directories(mesh INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(mesh PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

add_library(grid grid.f90 read.f90)
target_link_libraries(grid PUBLIC const reader mesh mpimod MPI::MPI_Fortran) # PUBLIC necessary for ifort, even with MPI target
target_include_directories(grid INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(grid PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
if(hdf5)
  target_sources(grid PRIVATE readgrid_hdf5.f90)
  target_link_libraries(grid PRIVATE h5fortran::h5fortran)
elseif(netcdf)
  target_sources(grid PRIVATE readgrid_nc4.f90)
  target_link_libraries(grid PRIVATE nc4fortran::nc4fortran)
else()
  target_sources(grid PRIVATE readgrid_raw.f90)
endif()
