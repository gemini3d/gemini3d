
add_executable(test_excess_mpi test_excess_mpi.f90)
target_link_libraries(test_excess_mpi PRIVATE autogrid)

add_test(NAME mpi:excessCPU COMMAND test_excess_mpi)
set_property(TEST mpi:excessCPU PROPERTY LABELS unit)
