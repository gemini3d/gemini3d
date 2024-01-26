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
if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  set_property(TARGET gemini_projection PROPERTY LINKER_LANGUAGE Fortran)
else()
  set_property(TARGET gemini_projection PROPERTY LINKER_LANGUAGE CXX)
endif()


add_test(NAME gemini_projection COMMAND gemini_projection)
set_tests_properties(gemini_projection PROPERTIES
LABELS unit
FIXTURES_REQUIRED mpi_fxt
)
test_mpi_launcher(gemini_projection 2)
