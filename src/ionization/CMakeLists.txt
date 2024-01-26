set_property(DIRECTORY PROPERTY  LABELS ionization)

add_library(fang OBJECT fang.f90)
target_link_libraries(fang PRIVATE const)

add_library(fang_run OBJECT fang_run.f90)
target_link_libraries(fang_run PRIVATE fang msis::msis const)

add_library(ionization OBJECT ionization.f90
$<TARGET_OBJECTS:grid>
)

if(glow)
  # must have trailing slash in path
  set(glow_data_dir "${glow_SOURCE_DIR}/data/")
  if(NOT EXISTS ${glow_data_dir}/ssflux_hint.dat)
    message(FATAL_ERROR "GLOW data directory is missing ssflux_hint.dat: ${glow_data_dir}")
  endif()

  configure_file(glow_run.in.f90 glow_run.f90 @ONLY)

  target_sources(ionization PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/glow_run.f90)
else()
  target_sources(ionization PRIVATE glow_dummy.f90)
endif()

target_link_libraries(ionization PRIVATE const meshobj timeutils fang
$<$<BOOL:${glow}>:glow::glow>
)

add_library(ionization_mpi OBJECT ionization_mpi.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(ionization_mpi PRIVATE ionization const neutral meshobj timeutils fang
mpimod
MPI::MPI_Fortran
)

add_library(precipBCs OBJECT boundary_conditions/precipBCs_mod.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(precipBCs PRIVATE reader const meshobj interp timeutils gemini3d_config
inputdata precipdata
MPI::MPI_Fortran
)
