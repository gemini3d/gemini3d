set_property(DIRECTORY PROPERTY LABEL neutral)

add_executable(test_proj test_proj.f90
)

foreach(t newton spherical geomagnetic meshobj meshobj_cart meshobj_dipole grid autogrid neutral interp
neutraldata neutraldata3D neutraldata3D_mpi neutraldata3D_geom_mpi neutraldata3D_geog_mpi neutraldata2D neutraldata2Daxisymm neutraldata2Dcart
)

  target_sources(test_proj PRIVATE $<TARGET_OBJECTS:${t}>)

endforeach()

target_include_directories(test_proj PRIVATE ${PROJECT_BINARY_DIR}/include)
target_link_libraries(test_proj PRIVATE
gemini3d_config const reader inputdata timeutils mpimod gemini3d_sysinfo
MPI::MPI_Fortran
msis::msis_ifc
$<$<BOOL:${hwm14}>:hwm14::hwm_ifc>
)
if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  set_property(TARGET test_proj PROPERTY LINKER_LANGUAGE Fortran)
else()
  set_property(TARGET test_proj PROPERTY LINKER_LANGUAGE CXX)
endif()


add_test(NAME Projection
COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_proj>
)
set_tests_properties(Projection PROPERTIES
PROCESSORS 2
LABELS unit
FIXTURES_REQUIRED mpi_fxt
)
if(DEFINED mpi_tmpdir)
  set_property(TEST Projection PROPERTY ENVIRONMENT "TMPDIR=${mpi_tmpdir}")
endif()

dll_test_path("ffilesystem::filesystem" Projection)
