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


# ---
if(NOT ${PROJECT_NAME}_BUILD_TESTING)
  return()
endif()

add_subdirectory(coord)
add_subdirectory(diffusion)
add_subdirectory(io)
add_subdirectory(ionization)
add_subdirectory(interpolation)
add_subdirectory(mpi)
add_subdirectory(neutral)
