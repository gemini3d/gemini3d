add_executable(test_formats test_formats.f90)
target_compile_options(test_formats PRIVATE ${FFLAGS})
target_link_libraries(test_formats PRIVATE timeutils)
set_target_properties(test_formats PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME DateFormats COMMAND test_formats)
set_tests_properties(DateFormats PROPERTIES TIMEOUT 5.0
  FIXTURES_SETUP IOfmt)

find_package(MPI REQUIRED COMPONENTS Fortran)

add_executable(testmpi test_mpi.f90)
target_link_libraries(testmpi PRIVATE MPI::MPI_Fortran)
target_compile_options(testmpi PRIVATE ${FFLAGS})
add_test(NAME MPIexist COMMAND testmpi)

set_tests_properties(MPIexist PROPERTIES
  TIMEOUT 5
  FIXTURES_SETUP MPIMUMPS)
#---------------


add_executable(testmumps test_mumps.f90)
target_include_directories(testmumps PRIVATE ${MUMPS_INCLUDE_DIRS})
target_compile_options(testmumps PRIVATE ${FFLAGS})
target_link_libraries(testmumps PRIVATE
  ${MUMPS_LIBRARIES} MPI::MPI_Fortran)

add_test(NAME MUMPS
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} $<TARGET_FILE:testmumps>
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

set_tests_properties(MUMPS PROPERTIES
  TIMEOUT 5
  DEPENDS MPIexist
  FIXTURES_SETUP MPIMUMPS)
