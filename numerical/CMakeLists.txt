add_library(const constants/phys_consts.F90)
target_compile_definitions(const PRIVATE REALBITS=${realbits})
target_compile_options(const PRIVATE ${FFLAGS})
set_target_properties(const PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

add_subdirectory(interpolation)

add_library(mpimod mpimod/mpimod.F90 mpimod/mpisend.f90 mpimod/mpirecv.f90 mpimod/mpihalo.f90)
target_link_libraries(mpimod PRIVATE const MPI::MPI_Fortran)
target_compile_definitions(mpimod PRIVATE REALBITS=${realbits})
target_compile_options(mpimod PRIVATE ${FFLAGS})
set_target_properties(mpimod PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

add_library(grid grid/grid.f90 grid/read.f90)
target_link_libraries(grid PUBLIC const mpimod MPI::MPI_Fortran) # PUBLIC necessary for ifort, even with MPI target
target_compile_options(grid PRIVATE ${FFLAGS})
set_target_properties(grid PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

add_library(calculus calculus/calculus.f90)
target_link_libraries(calculus PRIVATE grid)
target_compile_options(calculus PRIVATE ${FFLAGS})
set_target_properties(calculus PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

add_subdirectory(diffusion)

# -- more libs

add_library(advec advection/advec_mpi.f90)
target_link_libraries(advec PRIVATE const mpimod grid)
target_compile_options(advec PRIVATE ${FFLAGS})
set_target_properties(advec PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

add_subdirectory(potential)