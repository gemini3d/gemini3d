set_directory_properties(PROPERTIES LABELS ionization)

add_library(fang OBJECT fang.f90)
target_link_libraries(fang PRIVATE const)

add_library(fang_run OBJECT fang_run.f90)
target_link_libraries(fang_run PRIVATE fang msis_ifc const)

if(BUILD_TESTING)

add_executable(test_fang test_fang.f90)
target_link_libraries(test_fang PRIVATE fang fang_run assert const msis_ifc)
add_test(NAME unit:fangIonize COMMAND $<TARGET_FILE:test_fang>)
set_tests_properties(unit:fangIonize PROPERTIES
  TIMEOUT 10
  LABELS unit)

if(Python3_FOUND)
  add_test(NAME unit:fangIonize:python
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_fang.py $<TARGET_FILE:test_fang>)
  set_tests_properties(unit:fangIonize:python PROPERTIES
    TIMEOUT 30
    LABELS "unit;python")
endif(Python3_FOUND)

endif(BUILD_TESTING)

if(glow)
  include(${PROJECT_SOURCE_DIR}/cmake/ext_libs/glow.cmake)

  # DATADIR comes from glow project BUILD_INTERFACE COMPILE_DEFINITIONS
  # this does not seem to work with get_target_property()
  # so preprocess instead of configure_file()
  add_library(ionization OBJECT ionization.f90 glow_run.F90)
  target_link_libraries(ionization PUBLIC glow::glow)
else()
  add_library(ionization OBJECT ionization.f90 glow_dummy.f90)
endif(glow)
target_link_libraries(ionization PRIVATE const neutral grid meshobj timeutils fang mpimod MPI::MPI_Fortran)

add_library(precipBCs OBJECT boundary_conditions/precipBCs_mod.f90)
target_compile_options(precipBCs PRIVATE ${gcc10opts})
target_link_libraries(precipBCs PRIVATE mpimod reader const grid meshobj interp timeutils MPI::MPI_Fortran config)
