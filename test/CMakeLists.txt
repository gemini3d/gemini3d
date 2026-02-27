set_property(DIRECTORY PROPERTY LABELS standalone_libraries)

if(matlab)
  include(${PROJECT_SOURCE_DIR}/cmake/matlab.cmake)
endif(matlab)

add_compile_options(${gemini3d_flags})
# scope is this directory on down

# Libraries such as MPI, HDF5, Scalapack, MUMPS, are so important
# we test each standalone to mitigate confusing build or runtime errors
# if a particular library isn't working.

# --- HDF5
add_executable(test_hdf5_standalone_C test_hdf5.c)
target_link_libraries(test_hdf5_standalone_C PRIVATE h5fortran::h5fortran)
add_test(HDF5_standalone_C test_hdf5_standalone_C)
if(WIN32 AND CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  target_link_libraries(test_hdf5_standalone_C PRIVATE shlwapi)
endif()

add_executable(test_hdf5_standalone_Fortran test_hdf5.f90)
target_link_libraries(test_hdf5_standalone_Fortran PRIVATE h5fortran::h5fortran)
add_test(HDF5_standalone_Fortran test_hdf5_standalone_Fortran)
if(WIN32 AND CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  target_link_libraries(test_hdf5_standalone_Fortran PRIVATE shlwapi)
endif()

# --- MPI

add_executable(gemini_mpi_standalone test_mpi.f90)
target_link_libraries(gemini_mpi_standalone PRIVATE MPI::MPI_Fortran)

add_test(NAME GeminiMPIstandalone COMMAND gemini_mpi_standalone ${Ncpu})
test_mpi_launcher(gemini_mpi_standalone GeminiMPIstandalone ${Ncpu})

set_tests_properties(GeminiMPIstandalone PROPERTIES
FIXTURES_SETUP mpi_fxt
RESOURCE_LOCK cpu_mpi
)

# --- MUMPS
add_executable(gemini_mumps_standalone test_mumps.f90)
target_link_libraries(gemini_mumps_standalone PRIVATE mumps_interface
MUMPS::MUMPS MPI::MPI_Fortran
)
# same link logic applied to gemini.bin
target_include_directories(gemini_mumps_standalone PRIVATE ${MPI_Fortran_INCLUDE_DIRS})
# this include is for robustness so mpi_f08.mod isn't missed e.g. on Linux systems

add_test(NAME GeminiMUMPSstandalone COMMAND gemini_mumps_standalone)

set_tests_properties(GeminiMUMPSstandalone PROPERTIES
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
REQUIRED_FILES ${CMAKE_CURRENT_SOURCE_DIR}/input_simpletest_real.nml
FIXTURES_REQUIRED mpi_fxt
FIXTURES_SETUP mumps_fxt
RESOURCE_LOCK cpu_mpi
)
test_mpi_launcher(gemini_mumps_standalone GeminiMUMPSstandalone 2)



add_test(NAME gemini_fortran COMMAND gemini.bin -h)
set_tests_properties(gemini_fortran PROPERTIES
PASS_REGULAR_EXPRESSION "EOF: gemini.bin"
FIXTURES_SETUP gemini_exe_fxt
FIXTURES_REQUIRED "mumps_fxt"
LABELS unit
)


if(TARGET gemini_c.bin)

add_test(NAME gemini_cpp COMMAND gemini_c.bin -h)
set_tests_properties(gemini_cpp PROPERTIES
PASS_REGULAR_EXPRESSION "EOF: gemini.bin"
FIXTURES_SETUP gemini_exe_fxt
FIXTURES_REQUIRED "mumps_fxt"
LABELS "unit;Cpp"
)

endif()

add_test(NAME gemini3d_run_ok COMMAND gemini3d.run)
set_tests_properties(gemini3d_run_ok PROPERTIES
PASS_REGULAR_EXPRESSION "EOF: gemini3d.run"
FIXTURES_SETUP gemini_exe_fxt
)

add_subdirectory(compare)


if(NOT gemini3d_BUILD_TESTING)
  return()
endif()


add_executable(echo_path echo_path.f90)
add_test(NAME EchoEnvVar COMMAND echo_path)


add_test(NAME magcalc_exe_ok COMMAND magcalc.bin)
set_tests_properties(magcalc_exe_ok PROPERTIES
PASS_REGULAR_EXPRESSION "EOF: magcalc.bin"
LABELS unit
)

add_test(NAME magcalc_run_ok COMMAND magcalc.run)
set_tests_properties(magcalc_run_ok PROPERTIES
PASS_REGULAR_EXPRESSION "EOF: magcalc.run"
)


add_subdirectory(coord)
add_subdirectory(diffusion)
add_subdirectory(interpolation)
add_subdirectory(io)
add_subdirectory(ionization)
add_subdirectory(mpi)
add_subdirectory(neutral)
add_subdirectory(potential)
add_subdirectory(temporal)
