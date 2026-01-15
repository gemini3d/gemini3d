set_property(DIRECTORY PROPERTY  LABEL neutral)


add_library(neutral OBJECT neutral.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(neutral PRIVATE gemini3d_config const reader meshobj
inputdata
neutraldata
)


add_library(neutral_background OBJECT neutral_background_mod.f90 atmos.f90 wind.f90
$<$<NOT:$<BOOL:${hwm14}>>:${PROJECT_SOURCE_DIR}/src/vendor/hwm14_dummy.f90>
$<TARGET_OBJECTS:grid>
)
target_link_libraries(neutral_background PRIVATE gemini3d_config const reader meshobj interp mpimod grid
timeutils
inputdata
neutraldataBG
neutral
h5fortran::h5fortran
MPI::MPI_Fortran
msis::msis
$<$<BOOL:${hwm14}>:hwm14::hwm_ifc>
ffilesystem::filesystem
)


add_library(neutral_perturbations OBJECT neutral_perturbations_mod.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(neutral_perturbations PRIVATE gemini3d_config const reader meshobj interp mpimod grid
timeutils
inputdata
neutral neutraldata neutraldata3D neutraldata3D_fclaw neutraldata3D_mpi neutraldata3D_geog_mpi neutraldata3D_geom_mpi
neutraldata2D neutraldata2Dcart neutraldata2Daxisymm
h5fortran::h5fortran
MPI::MPI_Fortran
msis::msis
$<$<BOOL:${hwm14}>:hwm14::hwm_ifc>
ffilesystem::filesystem
)
