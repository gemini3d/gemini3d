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
target_include_directories(mpimod PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
# this include is for robustness so mpi_f08.mod isn't missed e.g. on Linux systems
target_link_libraries(mpimod PRIVATE const
MPI::MPI_Fortran
)
