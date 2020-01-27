add_library(const constants/phys_consts.F90)
target_include_directories(const INTERFACE ${CMAKE_CURRENT_BINARY_DIR})
target_compile_definitions(const PRIVATE REALBITS=${realbits})

add_subdirectory(interpolation)

add_library(mpimod mpimod/mpimod.F90 mpimod/mpisend.f90 mpimod/mpirecv.f90 mpimod/mpihalo.f90)
target_link_libraries(mpimod PRIVATE const MPI::MPI_Fortran)
target_compile_definitions(mpimod PRIVATE REALBITS=${realbits})

add_subdirectory(grid)

add_library(calculus calculus/calculus.f90 calculus/gradient.f90 calculus/integral.f90 calculus/div.f90)
target_link_libraries(calculus PRIVATE mesh)

add_subdirectory(diffusion)

# -- more libs

add_library(advec advection/advec_mpi.f90)
target_link_libraries(advec PRIVATE const mpimod grid)

add_subdirectory(potential)