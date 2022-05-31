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
