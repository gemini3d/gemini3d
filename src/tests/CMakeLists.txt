add_executable(test_mpi test_mpi.f90)
target_link_libraries(test_mpi MPI::MPI_Fortran)
add_test(NAME unit:mpi_basic COMMAND $<TARGET_FILE:test_mpi>)

set_tests_properties(unit:mpi_basic PROPERTIES
TIMEOUT 15
FIXTURES_SETUP MPIMUMPS)

# ---

if("d" IN_LIST arith)
  add_executable(gemini_scalapack_d test_scalapack_d.f90)
  target_link_libraries(gemini_scalapack_d ${SCALAPACK_LIBRARIES} ${LAPACK_LIBRARIES} MPI::MPI_Fortran)
  target_include_directories(gemini_scalapack_d PRIVATE ${SCALAPACK_INCLUDE_DIRS})

  add_test(NAME unit:scalapack_d COMMAND gemini_scalapack_d)
  set_tests_properties(unit:scalapack_d PROPERTIES
  TIMEOUT 15
  DEPENDS unit:mpi_basic
  FIXTURES_SETUP MPIMUMPS)
endif()

# ---
add_executable(test_mumps test_mumps.f90)
target_link_libraries(test_mumps mumps_interface MPI::MPI_Fortran)

add_test(NAME unit:mumps_basic
COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_mumps>
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

set_tests_properties(unit:mumps_basic PROPERTIES
TIMEOUT 15
DEPENDS unit:mpi_basic
PROCESSORS 2
RUN_SERIAL true
FIXTURES_SETUP MPIMUMPS)
