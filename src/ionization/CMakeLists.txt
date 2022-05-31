set_directory_properties(PROPERTIES LABELS ionization)

add_library(fang OBJECT fang.f90)
target_link_libraries(fang PRIVATE const)

add_library(fang_run OBJECT fang_run.f90)
target_link_libraries(fang_run PRIVATE fang msis::msis_ifc const)

if(BUILD_TESTING)

add_executable(test_fang test_fang.f90)
target_link_libraries(test_fang PRIVATE fang fang_run assert const msis::msis_ifc)

add_test(NAME fangIonize COMMAND test_fang)
set_tests_properties(fangIonize PROPERTIES
TIMEOUT 10
LABELS unit
)
dll_test_path(msis::msis_ifc fangIonize)

if(python)
  add_test(NAME fangIonize:python
  COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_fang.py $<TARGET_FILE:test_fang>
  )
  set_tests_properties(fangIonize:python PROPERTIES
  TIMEOUT 30
  LABELS "unit;python"
  )
endif()

endif(BUILD_TESTING)


add_library(ionization OBJECT ionization.f90)

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

target_link_libraries(ionization PRIVATE const grid meshobj timeutils fang
$<$<BOOL:${glow}>:glow::glow>
)

add_library(ionization_mpi OBJECT ionization_mpi.f90)
target_link_libraries(ionization_mpi PRIVATE ionization const neutral grid meshobj timeutils fang
mpimod
MPI::MPI_Fortran
)

add_library(precipBCs OBJECT boundary_conditions/precipBCs_mod.f90)
target_compile_options(precipBCs PRIVATE ${gfortran_opts})
target_link_libraries(precipBCs PRIVATE reader const grid meshobj interp timeutils config
inputdata precipdata
)
