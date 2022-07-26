set_property(DIRECTORY PROPERTY  LABELS unit)

if(matlab)
  include(${PROJECT_SOURCE_DIR}/cmake/matlab.cmake)
endif(matlab)

include(${PROJECT_SOURCE_DIR}/cmake/python.cmake)

if(mpi)
# keep this if(mpi) because these tests truly use MPI -- don't work with MPI stub

add_executable(test_mpi test_mpi.f90)
target_link_libraries(test_mpi PRIVATE MPI::MPI_Fortran)

add_test(NAME mpi_basic
COMMAND COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${Ncpu} $<TARGET_FILE:test_mpi> ${Ncpu}
)

set_tests_properties(mpi_basic PROPERTIES
TIMEOUT 30
PROCESSORS ${Ncpu}
FIXTURES_SETUP mpi_fxt
RESOURCE_LOCK cpu_mpi
)

# --- SCALAPACK

add_executable(gemini_scalapack test_scalapack.f90)
target_link_libraries(gemini_scalapack PRIVATE SCALAPACK::SCALAPACK LAPACK::LAPACK MPI::MPI_Fortran)

add_test(NAME scalapack COMMAND gemini_scalapack)
set_tests_properties(scalapack PROPERTIES
TIMEOUT 15
FIXTURES_REQUIRED mpi_fxt
FIXTURES_SETUP scalapack_fxt
)

# --- MUMPS
add_executable(test_mumps test_mumps.f90)
target_link_libraries(test_mumps PRIVATE mumps_interface
MUMPS::MUMPS
"$<$<BOOL:${MUMPS_HAVE_Scotch}>:Scotch::Scotch;METIS::METIS>"
)
# same link logic applied to gemini.bin

add_test(NAME mumps_basic
COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_mumps>
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

set_tests_properties(mumps_basic PROPERTIES
TIMEOUT 15
PROCESSORS 2
RESOURCE_LOCK cpu_mpi
REQUIRED_FILES ${CMAKE_CURRENT_SOURCE_DIR}/input_simpletest_real.nml
FIXTURES_REQUIRED scalapack_fxt
FIXTURES_SETUP mumps_fxt
)

dll_test_path("SCALAPACK::SCALAPACK;MUMPS::MUMPS" "scalapack;mumps_basic")


endif(mpi)


add_test(NAME gemini_fortran COMMAND gemini.bin -h)
set_tests_properties(gemini_fortran PROPERTIES
PASS_REGULAR_EXPRESSION "EOF: gemini.bin"
FIXTURES_SETUP gemini_exe_fxt
FIXTURES_REQUIRED "mumps_fxt"
LABELS unit
)
dll_test_path("ffilesystem::filesystem;h5fortran::h5fortran;LAPACK::LAPACK;HDF5::HDF5" gemini_fortran)


if(cpp)

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
TIMEOUT 10
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
TIMEOUT 10
)

add_test(NAME magcalc_run_ok COMMAND magcalc.run)
set_tests_properties(magcalc_run_ok PROPERTIES
PASS_REGULAR_EXPRESSION "EOF: magcalc.run"
TIMEOUT 10
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
