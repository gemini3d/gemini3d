add_library(temporal temporal.f90)
target_link_libraries(temporal PRIVATE const mesh mpimod MPI::MPI_Fortran)
target_compile_options(temporal PRIVATE ${FFLAGS})
set_target_properties(temporal PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

add_library(timeutils timeutils.f90)
target_link_libraries(timeutils PRIVATE const)
target_compile_options(timeutils PRIVATE ${FFLAGS})
set_target_properties(timeutils PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(test_formats test_formats.f90)
target_compile_options(test_formats PRIVATE ${FFLAGS})
target_link_libraries(test_formats PRIVATE timeutils const)
set_target_properties(test_formats PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME DateFormats COMMAND test_formats)
set_tests_properties(DateFormats PROPERTIES TIMEOUT 5.0
  FIXTURES_SETUP IOfmt)

add_executable(test_rollover test_dayrollover.f90)
target_compile_options(test_rollover PRIVATE ${FFLAGS})
target_link_libraries(test_rollover PRIVATE timeutils const)
set_target_properties(test_rollover PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME DateRollover COMMAND test_rollover)

add_executable(test_sza test_sza.f90)
target_compile_options(test_sza PRIVATE ${FFLAGS})
target_link_libraries(test_sza PRIVATE timeutils const)
set_target_properties(test_sza PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME SolarZenithAngle COMMAND test_sza)
