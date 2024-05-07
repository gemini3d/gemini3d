add_library(advec OBJECT advec.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(advec PRIVATE const meshobj)

add_library(advec_mpi OBJECT advec_mpi.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(advec_mpi PRIVATE const mpimod meshobj)

# apply compiler options
get_property(tgts DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)
foreach(t IN LISTS tgts)
  target_compile_options(${t} PRIVATE ${gemini3d_flags})
endforeach()
