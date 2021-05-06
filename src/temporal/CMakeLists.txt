add_library(timeutils OBJECT timeutils.f90)
target_link_libraries(timeutils PRIVATE const)

add_library(temporal OBJECT temporal.f90)
target_link_libraries(temporal PRIVATE config const meshobj mpimod MPI::MPI_Fortran)

# --- unit tests
if(${PROJECT_NAME}_BUILD_TESTING)

add_executable(test_formats test_formats.f90)
target_link_libraries(test_formats PRIVATE timeutils const)
add_test(NAME unit:DateFormats COMMAND $<TARGET_FILE:test_formats>)

add_executable(test_rollover test_dayrollover.f90)
target_link_libraries(test_rollover PRIVATE timeutils const)
add_test(NAME unit:DateRollover COMMAND $<TARGET_FILE:test_rollover>)

add_executable(test_sza test_sza.f90)
target_link_libraries(test_sza PRIVATE timeutils const)
add_test(NAME unit:SolarZenithAngle COMMAND $<TARGET_FILE:test_sza>)

endif()
