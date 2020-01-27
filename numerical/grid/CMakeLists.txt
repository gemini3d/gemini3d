add_library(grid grid.f90 read.f90)
target_link_libraries(grid PUBLIC const mesh mpimod MPI::MPI_Fortran) # PUBLIC necessary for ifort, even with MPI target
target_include_directories(grid INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(grid PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

add_library(mesh mesh.f90)
target_include_directories(mesh
  INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(mesh PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
target_link_libraries(mesh PRIVATE const)