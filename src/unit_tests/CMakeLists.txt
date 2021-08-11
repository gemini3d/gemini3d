set_directory_properties(PROPERTIES LABELS unit)

if(NOT BUILD_TESTING)
  return()
endif()

if(NOT mpi)
# keep this if(mpi) because these tests truly use MPI -- don't work with MPI stub
  return()
endif()

add_executable(test_mpi test_mpi.f90)
target_link_libraries(test_mpi PRIVATE MPI::MPI_Fortran)

add_test(NAME unit:mpi_basic
COMMAND COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${Ncpu} $<TARGET_FILE:test_mpi> ${Ncpu})

set_tests_properties(unit:mpi_basic PROPERTIES
TIMEOUT 15
PROCESSORS ${Ncpu}
DEPENDS unit:HWLOC
FIXTURES_SETUP mpi_fxt
RESOURCE_LOCK cpu_mpi)

# --- SCALAPACK

if("d" IN_LIST arith)
  set(lamch "pdlamch")
elseif("s" IN_LIST arith)
  set(lamch "pslamch")
endif()
configure_file(test_scalapack.in.f90 test_scalapack.f90 @ONLY)
add_executable(gemini_scalapack ${CMAKE_CURRENT_BINARY_DIR}/test_scalapack.f90)
target_link_libraries(gemini_scalapack PRIVATE SCALAPACK::SCALAPACK LAPACK::LAPACK MPI::MPI_Fortran)

add_test(NAME unit:scalapack COMMAND $<TARGET_FILE:gemini_scalapack>)
set_tests_properties(unit:scalapack PROPERTIES
TIMEOUT 15
FIXTURES_REQUIRED mpi_fxt
FIXTURES_SETUP scalapack_fxt)

# --- MUMPS
add_executable(test_mumps test_mumps.f90)
target_link_libraries(test_mumps PRIVATE mumps_interface
  MUMPS::MUMPS SCALAPACK::SCALAPACK LAPACK::LAPACK MPI::MPI_Fortran)
if(Scotch_FOUND)
  target_link_libraries(test_mumps PRIVATE Scotch::Scotch)
endif()
if(METIS_FOUND)
  target_link_libraries(test_mumps PRIVATE METIS::METIS)
endif()
# same link logic applied to gemini.bin

add_test(NAME unit:mumps_basic
COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_mumps>
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

set_tests_properties(unit:mumps_basic PROPERTIES
TIMEOUT 15
PROCESSORS 2
RESOURCE_LOCK cpu_mpi
REQUIRED_FILES ${CMAKE_CURRENT_SOURCE_DIR}/input_simpletest_real.nml
FIXTURES_REQUIRED scalapack_fxt
FIXTURES_SETUP mumps_fxt)
