set_property(DIRECTORY PROPERTY  LABELS inputdata)

add_library(inputdata OBJECT inputdataobj.f90)
target_link_libraries(inputdata PRIVATE const gemini3d_config meshobj meshobj_dipole interp timeutils)

add_library(precipdata OBJECT precipdataobj.f90)
target_link_libraries(precipdata PRIVATE const gemini3d_config meshobj interp timeutils inputdata
MPI::MPI_Fortran
)

add_library(efielddata OBJECT efielddataobj.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(efielddata PRIVATE const gemini3d_config meshobj interp timeutils inputdata)

add_library(neutraldata OBJECT neutraldataobj.f90)
target_link_libraries(neutraldata PRIVATE const inputdata)

add_library(neutraldata2D OBJECT neutraldata2Dobj.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(neutraldata2D PRIVATE
const gemini3d_config meshobj interp timeutils neutraldata inputdata reader
h5fortran::h5fortran
)

#add_library(neutraldata2Daxisymm_fclaw OBJECT neutraldata2Daxisymmobj_fclaw.f90
#$<TARGET_OBJECTS:grid>
#)
#target_link_libraries(neutraldata2Daxisymm_fclaw PRIVATE 
#const gemini3d_config meshobj interp timeutils neutraldata neutraldata2D inputdata
#)

add_library(neutraldata2Daxisymm OBJECT neutraldata2Daxisymmobj.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(neutraldata2Daxisymm PRIVATE
const gemini3d_config meshobj interp timeutils neutraldata neutraldata2D inputdata reader
h5fortran::h5fortran
)

add_library(neutraldata2Dcart OBJECT neutraldata2Dcartobj.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(neutraldata2Dcart PRIVATE
const gemini3d_config meshobj interp timeutils neutraldata neutraldata2D inputdata reader
h5fortran::h5fortran
)

add_library(neutraldata3D OBJECT neutraldata3Dobj.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(neutraldata3D PRIVATE 
const gemini3d_config meshobj interp timeutils neutraldata inputdata reader
)

add_library(neutraldata3D_fclaw OBJECT neutraldata3Dobj_fclaw.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(neutraldata3D_fclaw PRIVATE 
const gemini3d_config meshobj interp timeutils neutraldata neutraldata3D inputdata
)

add_library(neutraldata3D_fclaw_axisymm OBJECT neutraldata3Dobj_fclaw_axisymm.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(neutraldata3D_fclaw_axisymm PRIVATE 
const gemini3d_config meshobj interp timeutils neutraldata neutraldata3D neutraldata3D_fclaw inputdata
)

add_library(neutraldata3D_fclaw_3Dx OBJECT neutraldata3Dobj_fclaw_3Dx.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(neutraldata3D_fclaw_3Dx PRIVATE 
const gemini3d_config meshobj interp timeutils neutraldata neutraldata3D neutraldata3D_fclaw inputdata
)

add_library(neutraldata3D_mpi OBJECT neutraldata3Dobj_mpi.f90 neuslab_mpi.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(neutraldata3D_mpi PRIVATE
const gemini3d_config meshobj interp timeutils neutraldata neutraldata3D inputdata reader
h5fortran::h5fortran MPI::MPI_Fortran
)

add_library(neutraldata3D_geom_mpi OBJECT neutraldata3Dobj_geom_mpi.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(neutraldata3D_geom_mpi PRIVATE
const gemini3d_config meshobj interp timeutils neutraldata neutraldata3D neutraldata3D_mpi inputdata reader geomagnetic
h5fortran::h5fortran MPI::MPI_Fortran
)

add_library(neutraldata3D_geog_mpi OBJECT neutraldata3Dobj_geog_mpi.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(neutraldata3D_geog_mpi PRIVATE
const gemini3d_config meshobj interp timeutils neutraldata neutraldata3D neutraldata3D_mpi inputdata reader
h5fortran::h5fortran MPI::MPI_Fortran
)
