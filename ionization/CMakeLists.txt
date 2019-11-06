add_library(fang fang.f90)
target_link_libraries(fang PRIVATE const)
set_target_properties(fang PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

add_library(fang_run fang_run.f90)
target_link_libraries(fang_run PRIVATE fang msis const)
set_target_properties(fang_run PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(test_fang test_fang.f90)
target_link_libraries(test_fang PRIVATE fang_run const)
set_target_properties(test_fang PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

if(PythonOK)
add_test(NAME FangIonize
  COMMAND ${Python3_EXECUTABLE}  ${CMAKE_CURRENT_SOURCE_DIR}/test_fang.py $<TARGET_FILE:test_fang>)
set_tests_properties(FangIonize PROPERTIES
                     TIMEOUT 10)
endif(PythonOK)


add_library(ionization ionization.f90 glow_dummy.f90)
target_link_libraries(ionization PRIVATE const neutral grid timeutils fang mpimod MPI::MPI_Fortran)
set_target_properties(ionization PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

if(glow)
  include(${CMAKE_SOURCE_DIR}/cmake/glow.cmake)

  add_library(ionization_glow ionization.f90 glow_run.F90)
  target_link_libraries(ionization_glow PUBLIC cglow const neutral grid timeutils fang mpimod MPI::MPI_Fortran)
  set_target_properties(ionization_glow PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
  target_include_directories(ionization_glow PUBLIC ${GLOW_BINARY_DIR})
  target_link_directories(ionization_glow PUBLIC ${GLOW_BINARY_DIR})
  # double quotes needed in next line since it injects directly into Fortran code
  target_compile_definitions(ionization_glow PRIVATE DATADIR="${CMAKE_BINARY_DIR}/GLOW_proj-prefix/data/")

  # "ionization" is in this line to avoid intermittant quirks / race condition with ionization.mod
  # when ionization_glow and ionization may compile simultaneously.
  add_dependencies(ionization_glow GLOW_proj ionization)
endif(glow)

add_library(precipBCs boundary_conditions/precipBCs_mod.f90)
target_link_libraries(precipBCs PRIVATE mpimod const grid interp timeutils MPI::MPI_Fortran)
target_compile_options(precipBCs PRIVATE ${FFLAGS})
set_target_properties(precipBCs PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})