set_property(DIRECTORY PROPERTY  LABELS temporal)

add_library(timeutils OBJECT timeutils.f90)
target_link_libraries(timeutils PRIVATE const)

add_library(temporal_mpi OBJECT temporal_mpi.f90)
target_link_libraries(temporal_mpi PRIVATE gemini3d_config const meshobj mpimod
MPI::MPI_Fortran
)

add_library(temporal OBJECT temporal.f90)
target_link_libraries(temporal PRIVATE const)

