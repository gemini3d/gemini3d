add_executable(testmpi test_mpi.f90)
target_link_libraries(testmpi MPI::MPI_Fortran)
add_test(NAME mpi_basic COMMAND testmpi)

set_tests_properties(mpi_basic PROPERTIES
  TIMEOUT 15
  FIXTURES_SETUP MPIMUMPS)
# ---

add_executable(test_scalapack test_scalapack.f90)
target_link_libraries(test_scalapack ${SCALAPACK_LIBRARIES} ${LAPACK_LIBRARIES} MPI::MPI_Fortran)
target_include_directories(test_scalapack PRIVATE ${SCALAPACK_INCLUDE_DIRS})
add_test(NAME scalapack_basic COMMAND testmpi)

set_tests_properties(scalapack_basic PROPERTIES
  TIMEOUT 15
  FIXTURES_SETUP MPIMUMPS)

# ---
add_executable(test_mumps test_mumps.f90)
target_link_libraries(test_mumps mumps_interface MPI::MPI_Fortran)

add_test(NAME mumps_basic
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_mumps>
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

set_tests_properties(mumps_basic PROPERTIES
  TIMEOUT 15
  DEPENDS MPIexist
  RUN_SERIAL true
  FIXTURES_SETUP MPIMUMPS)
