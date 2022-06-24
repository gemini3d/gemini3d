set_directory_properties(PROPERTIES LABEL neutral)

add_executable(test_proj test_proj.f90
$<$<NOT:$<BOOL:${netcdf}>>:$<TARGET_OBJECTS:nc4fortran_dummy>>
)

foreach(t newton spherical geomagnetic meshobj meshobj_cart meshobj_dipole grid autogrid neutral interp
neutraldata neutraldata3D_mpi neutraldata3D_geom_mpi neutraldata3D_geog_mpi neutraldata2D neutraldata2Daxisymm neutraldata2Dcart
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
  set_target_properties(test_proj PROPERTIES LINKER_LANGUAGE Fortran)
else()
  set_target_properties(test_proj PROPERTIES LINKER_LANGUAGE CXX)
endif()

set(_cmd $<TARGET_FILE:test_proj>)
if(mpi)
  list(PREPEND _cmd ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2)
endif()

add_test(NAME Projection
COMMAND ${_cmd}
)
set_tests_properties(Projection PROPERTIES
TIMEOUT 90
LABELS unit
FIXTURES_REQUIRED mpi_fxt
)
