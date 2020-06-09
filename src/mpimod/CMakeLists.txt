if(NOT MPI_OK)
  return()
endif()

if(realbits EQUAL 32)
  set(mpi_realprec mpi_real)
else()
  set(mpi_realprec mpi_double_precision)
endif()

configure_file(mpimod.in.f90 mpimod.f90)

add_library(mpimod OBJECT
  ${CMAKE_CURRENT_BINARY_DIR}/mpimod.f90
  mpisend.f90 mpirecv.f90 mpihalo.f90)
target_link_libraries(mpimod PRIVATE const MPI::MPI_Fortran)

add_executable(test_excess_mpi test_excess_mpi.f90)
target_link_libraries(test_excess_mpi mpimod)
add_test(NAME unit:mpi:excessCPU COMMAND test_excess_mpi)
set_tests_properties(unit:mpi:excessCPU PROPERTIES TIMEOUT 10)
