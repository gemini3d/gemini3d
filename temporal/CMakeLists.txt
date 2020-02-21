add_library(temporal OBJECT temporal.f90)
target_link_libraries(temporal PRIVATE const mesh mpimod MPI::MPI_Fortran)

add_library(timeutils OBJECT timeutils.f90)
target_link_libraries(timeutils PRIVATE const)

# -- optional tests

add_executable(test_formats test_formats.f90)
target_link_libraries(test_formats PRIVATE timeutils const)
add_test(NAME unit:DateFormats COMMAND test_formats)

add_executable(test_rollover test_dayrollover.f90)
target_link_libraries(test_rollover PRIVATE timeutils const)
add_test(NAME unit:DateRollover COMMAND test_rollover)

add_executable(test_sza test_sza.f90)
target_link_libraries(test_sza PRIVATE timeutils const)
add_test(NAME unit:SolarZenithAngle COMMAND test_sza)
