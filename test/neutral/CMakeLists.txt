set_property(DIRECTORY PROPERTY LABEL neutral)

add_executable(gemini_projection test_proj.f90
)

foreach(t IN ITEMS newton spherical geomagnetic meshobj meshobj_cart meshobj_dipole grid autogrid neutral interp
neutraldata neutraldata3D neutraldata3D_mpi neutraldata3D_geom_mpi neutraldata3D_geog_mpi neutraldata2D neutraldata2Daxisymm neutraldata2Dcart
)

  target_sources(gemini_projection PRIVATE $<TARGET_OBJECTS:${t}>)

endforeach()

target_include_directories(gemini_projection PRIVATE ${PROJECT_BINARY_DIR}/include)
target_link_libraries(gemini_projection PRIVATE
gemini3d_config const reader inputdata timeutils mpimod gemini3d_sysinfo
MPI::MPI_Fortran
msis::msis
$<$<BOOL:${hwm14}>:hwm14::hwm_ifc>
)
set_property(TARGET gemini_projection PROPERTY LINKER_LANGUAGE Fortran)
# for oneAPI

add_test(NAME GeminiProjection COMMAND gemini_projection)
set_tests_properties(GeminiProjection PROPERTIES
LABELS unit
FIXTURES_REQUIRED mpi_fxt
)
test_mpi_launcher(gemini_projection GeminiProjection 2)
