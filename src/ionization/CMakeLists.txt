set_property(DIRECTORY PROPERTY  LABELS ionization)

add_library(fang OBJECT fang.f90)
target_link_libraries(fang PRIVATE const)

add_library(fang_run OBJECT fang_run.f90)
target_link_libraries(fang_run PRIVATE fang msis::msis_ifc const)

add_library(ionization OBJECT ionization.f90
$<TARGET_OBJECTS:grid>
)

if(glow)
  get_target_property(_def glow::glow INTERFACE_COMPILE_DEFINITIONS)
  string(REGEX MATCH "^DATADIR\=\"(.+)\"$" _def "${_def}")
  set(glow_data_dir "${CMAKE_MATCH_1}")
  if(NOT IS_DIRECTORY "${glow_data_dir}")
    message(FATAL_ERROR "GLOW data directory not found: ${_def}
    ${glow_data_dir}")
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
