set_property(DIRECTORY PROPERTY LABELS standalone_libraries)

if(matlab)
  include(${PROJECT_SOURCE_DIR}/cmake/matlab.cmake)
endif(matlab)

# Libraries such as MPI, HDF5, Scalapack, MUMPS, are so important
# we test each standalone to mitigate confusing build or runtime errors
# if a particular library isn't working.

# --- HDF5
add_executable(test_hdf5_standalone_C test_hdf5.c)
target_link_libraries(test_hdf5_standalone_C PRIVATE HDF5::HDF5)
add_test(HDF5_standalone_C test_hdf5_standalone_C)

add_executable(test_hdf5_standalone_Fortran test_hdf5.f90)
target_link_libraries(test_hdf5_standalone_Fortran PRIVATE HDF5::HDF5)
add_test(HDF5_standalone_Fortran test_hdf5_standalone_Fortran)

# --- MPI

add_executable(test_mpi_standalone test_mpi.f90)
target_link_libraries(test_mpi_standalone PRIVATE MPI::MPI_Fortran)

add_test(NAME MPI_standalone
COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${Ncpu} $<TARGET_FILE:test_mpi_standalone> ${Ncpu}
)

set_tests_properties(MPI_standalone PROPERTIES
PROCESSORS ${Ncpu}
FIXTURES_SETUP mpi_fxt
)

# --- SCALAPACK

add_executable(test_scalapack_gemini test_scalapack.f90)
target_link_libraries(test_scalapack_gemini PRIVATE SCALAPACK::SCALAPACK LAPACK::LAPACK MPI::MPI_Fortran)

add_test(NAME SCALAPACK_standalone COMMAND test_scalapack_gemini)
set_tests_properties(SCALAPACK_standalone PROPERTIES
FIXTURES_REQUIRED mpi_fxt
FIXTURES_SETUP scalapack_fxt
)

# --- MUMPS
add_executable(MUMPS_standalone test_mumps.f90)
target_link_libraries(MUMPS_standalone PRIVATE mumps_interface
MUMPS::MUMPS
)
# same link logic applied to gemini.bin
target_include_directories(MUMPS_standalone PRIVATE ${MPI_Fortran_INCLUDE_DIRS})
# this include is for robustness so mpi_f08.mod isn't missed e.g. on Linux systems

add_test(NAME MUMPS_standalone
COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:MUMPS_standalone>
)

set_tests_properties(MUMPS_standalone PROPERTIES
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
PROCESSORS 2
REQUIRED_FILES ${CMAKE_CURRENT_SOURCE_DIR}/input_simpletest_real.nml
FIXTURES_REQUIRED scalapack_fxt
FIXTURES_SETUP mumps_fxt
)

set_property(TEST MPI_standalone MUMPS_standalone PROPERTY RESOURCE_LOCK cpu_mpi)
if(DEFINED mpi_tmpdir)
  set_property(TEST MPI_standalone MUMPS_standalone PROPERTY ENVIRONMENT "TMPDIR=${mpi_tmpdir}")
endif()

dll_test_path("SCALAPACK::SCALAPACK;MUMPS::MUMPS" "SCALAPACK_standalone;MUMPS_standalone")


add_test(NAME gemini_fortran COMMAND gemini.bin -h)
set_tests_properties(gemini_fortran PROPERTIES
PASS_REGULAR_EXPRESSION "EOF: gemini.bin"
FIXTURES_SETUP gemini_exe_fxt
FIXTURES_REQUIRED "mumps_fxt"
LABELS unit
)
dll_test_path("ffilesystem::filesystem;h5fortran::h5fortran;LAPACK::LAPACK;HDF5::HDF5" gemini_fortran)


if(TARGET gemini_c.bin)

add_test(NAME gemini_cpp COMMAND gemini_c.bin -h)
set_tests_properties(gemini_cpp PROPERTIES
PASS_REGULAR_EXPRESSION "EOF: gemini.bin"
FIXTURES_SETUP gemini_exe_fxt
FIXTURES_REQUIRED "mumps_fxt"
LABELS unit
)
dll_test_path("ffilesystem::filesystem;h5fortran::h5fortran;LAPACK::LAPACK;HDF5::HDF5" gemini_cpp)

endif()

add_test(NAME gemini3d_run_ok COMMAND gemini3d.run)
set_tests_properties(gemini3d_run_ok PROPERTIES
PASS_REGULAR_EXPRESSION "EOF: gemini3d.run"
FIXTURES_SETUP gemini_exe_fxt
)

# --- Windows needs DLLs on path
dll_test_path("ffilesystem::filesystem;h5fortran::h5fortran;HDF5::HDF5" "gemini3d_run_ok")


add_subdirectory(compare)


if(NOT ${PROJECT_NAME}_BUILD_TESTING)
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

dll_test_path("ffilesystem::filesystem;h5fortran::h5fortran;HDF5::HDF5" "magcalc_exe_ok;magcalc_run_ok")



add_subdirectory(coord)
add_subdirectory(diffusion)
add_subdirectory(interpolation)
add_subdirectory(io)
add_subdirectory(ionization)
add_subdirectory(mpi)
add_subdirectory(neutral)
add_subdirectory(potential)
add_subdirectory(temporal)
