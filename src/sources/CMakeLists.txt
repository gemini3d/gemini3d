add_library(sources OBJECT sources.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(sources PRIVATE calculus collision const meshobj grid)
add_library(sources_mpi OBJECT sources_mpi.f90)
target_link_libraries(sources_mpi PRIVATE mpimod)

# apply compiler options
get_property(tgts DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)
foreach(t IN LISTS tgts)
  target_compile_options(${t} PRIVATE ${gemini3d_flags})
endforeach()
