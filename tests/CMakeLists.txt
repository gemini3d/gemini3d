find_package(MPI REQUIRED COMPONENTS Fortran)

add_executable(testmpi test_mpi.f90)
target_link_libraries(testmpi PRIVATE MPI::MPI_Fortran)
target_compile_options(testmpi PRIVATE ${FFLAGS})
add_test(NAME MPIbasic COMMAND testmpi)

set_tests_properties(MPIbasic PROPERTIES
  TIMEOUT 15
  FIXTURES_SETUP MPIMUMPS)
# ---

add_executable(test_scalapack test_scalapack.f90)
target_link_libraries(test_scalapack PRIVATE ${SCALAPACK_LIBRARIES} ${LAPACK_LIBRARIES} MPI::MPI_Fortran)
target_include_directories(test_scalapack PRIVATE ${SCALAPACK_INCLUDE_DIRS})
target_compile_options(test_scalapack PRIVATE ${FFLAGS})
add_test(NAME ScalapackBasic COMMAND testmpi)

set_tests_properties(ScalapackBasic PROPERTIES
  TIMEOUT 15
  FIXTURES_SETUP MPIMUMPS)

# ---
add_executable(test_mumps test_mumps.F90)
target_include_directories(test_mumps PRIVATE ${MUMPS_INCLUDE_DIRS})
target_compile_options(test_mumps PRIVATE ${FFLAGS})
target_compile_definitions(test_mumps PRIVATE REALBITS=${realbits})
target_link_libraries(test_mumps PRIVATE ${MUMPS_LIBRARIES} MPI::MPI_Fortran)

if(mumps_external)
  add_dependencies(test_mumps MUMPS)
endif()

add_test(NAME MUMPSbasic
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_mumps>
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

set_tests_properties(MUMPSbasic PROPERTIES
  TIMEOUT 15
  DEPENDS MPIexist
  RUN_SERIAL true
  FIXTURES_SETUP MPIMUMPS)
