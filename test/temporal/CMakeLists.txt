set_property(DIRECTORY PROPERTY LABELS temporal)

add_executable(test_formats test_formats.f90)
target_link_libraries(test_formats PRIVATE timeutils const)
add_test(NAME DateFormats COMMAND test_formats)

add_executable(test_rollover test_dayrollover.f90)
target_link_libraries(test_rollover PRIVATE timeutils const)
add_test(NAME DateRollover COMMAND test_rollover)

add_executable(test_sza test_sza.f90)
target_link_libraries(test_sza PRIVATE timeutils const)
add_test(NAME SolarZenithAngle COMMAND test_sza)

set_property(TEST DateFormats DateRollover SolarZenithAngle PROPERTY LABELS unit)
