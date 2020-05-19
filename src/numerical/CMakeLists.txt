# --- generate constants
configure_file(constants/phys_consts.in.f90 phys_consts.f90)
add_library(const OBJECT ${CMAKE_CURRENT_BINARY_DIR}/phys_consts.f90)

# --- interpolation
add_subdirectory(interpolation)

# --- grid
add_subdirectory(grid)

# --- calculus
add_library(calculus OBJECT
  calculus/calculus.f90 calculus/gradient.f90 calculus/integral.f90 calculus/div.f90)
target_link_libraries(calculus PRIVATE mesh)

if(NOT MPI_OK)
  set_targ_props(calculus const interp mesh)
  return()
endif()

# --- mpi
if(realbits EQUAL 32)
  set(mpi_realprec mpi_real)
else()
  set(mpi_realprec mpi_double_precision)
endif()

configure_file(mpimod/mpimod.in.f90 mpimod.f90)

add_library(mpimod OBJECT
  ${CMAKE_CURRENT_BINARY_DIR}/mpimod.f90
  mpimod/mpisend.f90 mpimod/mpirecv.f90 mpimod/mpihalo.f90)
target_link_libraries(mpimod PRIVATE const MPI::MPI_Fortran)

add_subdirectory(diffusion)

# -- more libs

add_library(advec OBJECT advection/advec_mpi.f90)
target_link_libraries(advec PRIVATE const mpimod grid mesh)

add_subdirectory(potential)

set_targ_props(advec calculus const diffusion grid interp mesh mpimod mumps_interface PDEelliptic PDEparabolic potential)