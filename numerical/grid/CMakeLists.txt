add_library(grid grid.f90 mesh.f90 read.f90)
target_link_libraries(grid PUBLIC const mesh mpimod MPI::MPI_Fortran) # PUBLIC necessary for ifort, even with MPI target
target_compile_options(grid PRIVATE ${FFLAGS})
set_target_properties(grid PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

add_library(mesh mesh.f90)
target_link_libraries(mesh PRIVATE const)
target_compile_options(mesh PRIVATE ${FFLAGS})
set_target_properties(mesh PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})