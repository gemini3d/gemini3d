add_library(neutral OBJECT neutral.f90 atmos.f90 wind.f90
$<$<NOT:$<BOOL:${hwm14}>>:${PROJECT_SOURCE_DIR}/src/vendor/hwm14_dummy.f90>
)
target_compile_options(neutral PRIVATE ${gfortran_opts})
target_link_libraries(neutral PRIVATE config const reader grid meshobj interp mpimod
timeutils
inputdata neutraldata neutraldata3D_mpi
neutraldata2D neutraldata2Dcart neutraldata2Daxisymm
h5fortran::h5fortran
MPI::MPI_Fortran
msis::msis_ifc
$<$<BOOL:${hwm14}>:hwm14::hwm_ifc>
ffilesystem::filesystem
)

add_library(neutral_perturbations OBJECT neutral_perturbations.f90)
target_compile_options(neutral_perturbations PRIVATE ${gfortran_opts})
target_link_libraries(neutral_perturbations PRIVATE config const reader grid meshobj interp mpimod
timeutils
inputdata
neutral neutraldata neutraldata3D_mpi neutraldata3D_geog_mpi neutraldata3D_geom_mpi
neutraldata2D neutraldata2Dcart neutraldata2Daxisymm
h5fortran::h5fortran
MPI::MPI_Fortran
msis::msis_ifc
$<$<BOOL:${hwm14}>:hwm14::hwm_ifc>
ffilesystem::filesystem
)

# self-tests
if(${PROJECT_NAME}_BUILD_TESTING)

add_executable(test_proj test_proj.f90
$<$<NOT:$<BOOL:${netcdf}>>:$<TARGET_OBJECTS:nc4fortran_dummy>>
)

foreach(t newton spherical geomagnetic meshobj meshobj_cart meshobj_dipole grid autogrid neutral interp
neutraldata neutraldata3D_mpi neutraldata3D_geom_mpi neutraldata3D_geog_mpi neutraldata2D neutraldata2Daxisymm neutraldata2Dcart
)

  target_sources(test_proj PRIVATE $<TARGET_OBJECTS:${t}>)

endforeach()

target_include_directories(test_proj PRIVATE ${PROJECT_BINARY_DIR}/include)
target_link_libraries(test_proj PRIVATE config const reader inputdata
timeutils mpimod const
MPI::MPI_Fortran
msis::msis_ifc
$<$<BOOL:${hwm14}>:hwm14::hwm_ifc>
ffilesystem::filesystem
)
if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  set_target_properties(test_proj PROPERTIES LINKER_LANGUAGE Fortran)
else()
  set_target_properties(test_proj PROPERTIES LINKER_LANGUAGE CXX)
endif()

endif(${PROJECT_NAME}_BUILD_TESTING)
