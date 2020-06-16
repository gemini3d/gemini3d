# --- generate constants
configure_file(constants/phys_consts.in.f90 phys_consts.f90)
add_library(const OBJECT ${CMAKE_CURRENT_BINARY_DIR}/phys_consts.f90)

# --- interpolation
add_subdirectory(interpolation)

# --- grid
add_subdirectory(grid)

# --- calculus
configure_file(calculus/calculus.in.f90 calculus.f90)
add_library(calculus OBJECT
  ${CMAKE_CURRENT_BINARY_DIR}/calculus.f90 calculus/gradient.f90 calculus/integral.f90 calculus/div.f90)
target_link_libraries(calculus PRIVATE mesh)

if(NOT MPI_OK)
  set_targ_props(calculus const interp mesh)
  return()
endif()

add_subdirectory(diffusion)

# -- more libs

add_library(advec OBJECT advection/advec_mpi.f90)
target_link_libraries(advec PRIVATE const mpimod grid mesh)

add_subdirectory(potential)

set_targ_props(advec calculus const diffusion grid interp mesh mumps_interface PDEelliptic PDEparabolic potential)
