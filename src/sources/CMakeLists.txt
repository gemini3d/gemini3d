add_library(sources OBJECT sources.f90
$<TARGET_OBJECTS:grid> $<TARGET_OBJECTS:gemini3d_config>
)
target_link_libraries(sources PRIVATE calculus collision const meshobj grid)

# A callable conservative pair operator for qualification and later coupling.
# It is not substituted into the existing split source driver implicitly.
add_library(pair_exchange OBJECT pair_exchange.f90)
target_link_libraries(pair_exchange PRIVATE const)

add_library(sources_mpi OBJECT sources_mpi.f90)
target_link_libraries(sources_mpi PRIVATE mpimod)
