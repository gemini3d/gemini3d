add_library(sources OBJECT sources.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(sources PRIVATE calculus collision const meshobj grid)
add_library(sources_mpi OBJECT sources_mpi.f90)
target_link_libraries(sources_mpi PRIVATE mpimod)
