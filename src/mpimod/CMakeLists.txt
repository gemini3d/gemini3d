add_library(autogrid autogrid.f90)

if(realbits EQUAL 32)
  set(mpi_realprec mpi_real)
else()
  set(mpi_realprec mpi_double_precision)
endif()

configure_file(mpimod.in.f90 mpimod.f90 @ONLY)

add_library(mpimod OBJECT
  ${CMAKE_CURRENT_BINARY_DIR}/mpimod.f90
  mpisend.f90 mpirecv.f90 mpihalo.f90)
target_link_libraries(mpimod PRIVATE autogrid const MPI::MPI_Fortran)
if(NOT mpi)
  # the reason for this is the MPI-1 Fortran 77 linking in MUMPS
  # this seemed to be an easy way to avoid symbol clashes
  target_link_libraries(mpimod PRIVATE MUMPS::mpiseq)
endif()

add_executable(mpi_cli get_mpi_count.f90)
target_link_libraries(mpi_cli PRIVATE autogrid reader const pathlib timeutils config)
set_target_properties(mpi_cli PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

if(${PROJECT_NAME}_BUILD_TESTING)

add_executable(test_excess_mpi test_excess_mpi.f90)
target_link_libraries(test_excess_mpi PRIVATE autogrid)
add_test(NAME unit:mpi:excessCPU COMMAND $<TARGET_FILE:test_excess_mpi>)
set_tests_properties(unit:mpi:excessCPU PROPERTIES TIMEOUT 10)

endif()
