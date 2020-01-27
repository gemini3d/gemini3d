add_library(grid grid.f90 mesh.f90 read.f90)
target_link_libraries(grid PUBLIC const mesh mpimod MPI::MPI_Fortran) # PUBLIC necessary for ifort, even with MPI target

add_library(mesh mesh.f90)
target_include_directories(mesh INTERFACE ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(mesh PRIVATE const)