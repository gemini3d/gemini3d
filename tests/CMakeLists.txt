find_package(MPI REQUIRED COMPONENTS Fortran)

add_executable(testmpi test_mpi.f90)
target_link_libraries(testmpi PRIVATE MPI::MPI_Fortran)
target_compile_options(testmpi PRIVATE ${FFLAGS})
add_test(NAME MPIexist COMMAND testmpi)

set_tests_properties(MPIexist PROPERTIES
  TIMEOUT 5
  FIXTURES_SETUP MPIMUMPS)
#---------------


add_executable(testmumps test_mumps.F90)
target_include_directories(testmumps PRIVATE ${MUMPS_INCLUDE_DIRS})
target_compile_options(testmumps PRIVATE ${FFLAGS})
target_compile_definitions(testmumps PRIVATE REALBITS=${realbits})
target_link_libraries(testmumps PRIVATE
  ${MUMPS_LIBRARIES} MPI::MPI_Fortran)

add_test(NAME MUMPS
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} $<TARGET_FILE:testmumps>
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

set_tests_properties(MUMPS PROPERTIES
  TIMEOUT 5
  DEPENDS MPIexist
  FIXTURES_SETUP MPIMUMPS)
