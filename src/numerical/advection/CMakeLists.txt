add_library(transport_audit OBJECT transport_audit.f90)
target_link_libraries(transport_audit PRIVATE const meshobj)
add_library(advec OBJECT advec.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(advec PRIVATE const meshobj transport_audit)

add_library(advec_mpi OBJECT advec_mpi.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(advec_mpi PRIVATE const mpimod meshobj)
