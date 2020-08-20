add_executable(test_mpi test_mpi.f90)
target_link_libraries(test_mpi MPI::MPI_Fortran)
add_test(NAME unit:mpi_basic
COMMAND COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} $<TARGET_FILE:test_mpi> ${MPIEXEC_MAX_NUMPROCS})

set_tests_properties(unit:mpi_basic PROPERTIES
TIMEOUT 15
FIXTURES_SETUP mpi_fixture)

if(no_mpi)
return()
endif()

# --- SCALAPACK

if("d" IN_LIST arith)
  set(lamch "pdlamch")
elseif("s" IN_LIST arith)
  set(lamch "pslamch")
endif()
configure_file(test_scalapack.in.f90 test_scalapack.f90)
add_executable(gemini_scalapack ${CMAKE_CURRENT_BINARY_DIR}/test_scalapack.f90)
target_link_libraries(gemini_scalapack SCALAPACK::SCALAPACK LAPACK::LAPACK MPI::MPI_Fortran)

add_test(NAME unit:scalapack COMMAND gemini_scalapack)
set_tests_properties(unit:scalapack PROPERTIES
TIMEOUT 15
FIXTURES_REQUIRED mpi_fixture
FIXTURES_SETUP scalapack_fixture)

# --- MUMPS
add_executable(test_mumps test_mumps.f90)
target_link_libraries(test_mumps mumps_interface SCALAPACK::SCALAPACK LAPACK::LAPACK MPI::MPI_Fortran)

add_test(NAME unit:mumps_basic
COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_mumps>
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

set_tests_properties(unit:mumps_basic PROPERTIES
TIMEOUT 15
PROCESSORS 2
SKIP_RETURN_CODE 77
RUN_SERIAL true
FIXTURES_REQUIRED scalapack_fixture
FIXTURES_SETUP mumps_fixture)
