
add_executable(test_excess_mpi test_excess_mpi.f90)
target_link_libraries(test_excess_mpi PRIVATE autogrid)

add_test(NAME mpi:excessCPU COMMAND test_excess_mpi)
set_property(TEST mpi:excessCPU PROPERTY LABELS unit)

# apply compiler options
get_property(tgts DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)
foreach(t IN LISTS tgts)
  target_compile_options(${t} PRIVATE ${gemini3d_flags})
endforeach()
