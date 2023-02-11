set_property(DIRECTORY PROPERTY  LABEL neutral)

add_library(neutral OBJECT neutral.f90 atmos.f90 wind.f90
$<$<NOT:$<BOOL:${hwm14}>>:${PROJECT_SOURCE_DIR}/src/vendor/hwm14_dummy.f90>
$<TARGET_OBJECTS:grid>
)
target_link_libraries(neutral PRIVATE gemini3d_config const reader meshobj interp mpimod
timeutils
inputdata neutraldata neutraldata3D_mpi
neutraldata2D neutraldata2Dcart neutraldata2Daxisymm
h5fortran::h5fortran
MPI::MPI_Fortran
msis::msis_ifc
$<$<BOOL:${hwm14}>:hwm14::hwm_ifc>
ffilesystem::filesystem
)

add_library(neutral_perturbations OBJECT neutral_perturbations.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(neutral_perturbations PRIVATE gemini3d_config const reader meshobj interp mpimod
timeutils
inputdata
neutral neutraldata neutraldata3D neutraldata3D_fclaw neutraldata3D_mpi neutraldata3D_geog_mpi neutraldata3D_geom_mpi
neutraldata2D neutraldata2Dcart neutraldata2Daxisymm
h5fortran::h5fortran
MPI::MPI_Fortran
msis::msis_ifc
$<$<BOOL:${hwm14}>:hwm14::hwm_ifc>
ffilesystem::filesystem
)
