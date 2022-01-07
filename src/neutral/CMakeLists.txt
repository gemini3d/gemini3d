add_library(neutral OBJECT neutral.f90 atmos.f90 wind.f90)
target_compile_options(neutral PRIVATE ${gfortran_opts})
target_link_libraries(neutral PRIVATE config const reader grid meshobj interp mpimod
timeutils
inputdata neutraldata neutraldata3D_mpi
neutraldata2D neutraldata2Dcart neutraldata2Daxisymm
h5fortran::h5fortran
MPI::MPI_Fortran
msis::msis_ifc
hwm14::hwm_ifc
ffilesystem::filesystem ${lib_filesystem}
)

# self-tests
add_executable(test_proj test_proj.f90)

foreach(t newton spherical geomagnetic meshobj meshobj_cart meshobj_dipole grid autogrid neutral interp
neutraldata neutraldata3D neutraldata2D neutraldata2Dcart neutraldata2Daxisymm
)

  target_sources(test_proj PRIVATE $<TARGET_OBJECTS:${t}>)

endforeach()

target_include_directories(test_proj PRIVATE ${PROJECT_BINARY_DIR}/include)
target_link_libraries(test_proj PRIVATE config const reader inputdata
timeutils mpimod const
MPI::MPI_Fortran
msis::msis_ifc
hwm14::hwm_ifc
ffilesystem::filesystem ${lib_filesystem}
)
if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  set_target_properties(test_proj PROPERTIES LINKER_LANGUAGE Fortran)
else()
  set_target_properties(test_proj PROPERTIES LINKER_LANGUAGE CXX)
endif()
