add_library(autogrid OBJECT autogrid.f90)

if(realbits EQUAL 32)
  set(mpi_realprec mpi_real)
else()
  set(mpi_realprec mpi_double_precision)
endif()

configure_file(mpimod.in.f90 mpimod.f90 @ONLY)

add_library(mpimod OBJECT
${CMAKE_CURRENT_BINARY_DIR}/mpimod.f90
$<TARGET_OBJECTS:autogrid>
mpisend.f90 mpirecv.f90 mpihalo.f90
)
target_link_libraries(mpimod PRIVATE const
MPI::MPI_Fortran
$<$<NOT:$<BOOL:${mpi}>>:MUMPS::MPISEQ>
)
# MPI-1 Fortran 77 linking in MUMPS--easy way to avoid symbol clashes

#-----

if(BUILD_TESTING)

add_executable(test_excess_mpi test_excess_mpi.f90)
target_link_libraries(test_excess_mpi PRIVATE autogrid)

add_test(NAME mpi:excessCPU COMMAND test_excess_mpi)
set_tests_properties(mpi:excessCPU PROPERTIES
TIMEOUT 10
LABELS unit
)
endif()
