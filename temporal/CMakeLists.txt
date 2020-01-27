add_library(temporal temporal.f90)
target_link_libraries(temporal PRIVATE const mesh mpimod MPI::MPI_Fortran)
target_include_directories(temporal INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(temporal PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

add_library(timeutils timeutils.f90)
target_link_libraries(timeutils PRIVATE const)
target_include_directories(timeutils INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(timeutils PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

# -- optional tests

add_executable(test_formats test_formats.f90)
target_link_libraries(test_formats PRIVATE timeutils const)
add_test(NAME DateFormats COMMAND test_formats)
set_tests_properties(DateFormats PROPERTIES TIMEOUT 5.0
  FIXTURES_SETUP IOfmt)

add_executable(test_rollover test_dayrollover.f90)
target_link_libraries(test_rollover PRIVATE timeutils const)
add_test(NAME DateRollover COMMAND test_rollover)

add_executable(test_sza test_sza.f90)
target_link_libraries(test_sza PRIVATE timeutils const)
add_test(NAME SolarZenithAngle COMMAND test_sza)
