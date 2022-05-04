set_directory_properties(PROPERTIES LABELS inputdata)

add_library(inputdata OBJECT inputdataobj.f90)
target_link_libraries(inputdata PRIVATE const config meshobj meshobj_dipole interp timeutils)

add_library(precipdata OBJECT precipdataobj.f90)
target_link_libraries(precipdata PRIVATE const config meshobj interp timeutils inputdata MPI::MPI_Fortran)

add_library(efielddata OBJECT efielddataobj.f90)
target_link_libraries(efielddata PRIVATE const config grid meshobj interp timeutils inputdata)

add_library(neutraldata OBJECT neutraldataobj.f90)
target_link_libraries(neutraldata PRIVATE const inputdata)

add_library(neutraldata2D OBJECT neutraldata2Dobj.f90)
target_link_libraries(neutraldata2D PRIVATE const config grid meshobj interp timeutils neutraldata inputdata reader h5fortran::h5fortran)

add_library(neutraldata2Daxisymm OBJECT neutraldata2Daxisymmobj.f90)
target_link_libraries(neutraldata2Daxisymm PRIVATE const config grid meshobj interp timeutils neutraldata neutraldata2D inputdata reader h5fortran::h5fortran)
target_compile_options(neutraldata2Daxisymm PRIVATE ${gfortran_opts})

add_library(neutraldata2Dcart OBJECT neutraldata2Dcartobj.f90)
target_link_libraries(neutraldata2Dcart PRIVATE const config grid meshobj interp timeutils neutraldata neutraldata2D inputdata reader h5fortran::h5fortran)
target_compile_options(neutraldata2Dcart PRIVATE ${gfortran_opts})

add_library(neutraldata3D_mpi OBJECT neutraldata3Dobj_mpi.f90 neuslab_mpi.f90)
target_link_libraries(neutraldata3D_mpi PRIVATE const config grid meshobj interp timeutils neutraldata inputdata reader MPI::MPI_Fortran h5fortran::h5fortran)
target_compile_options(neutraldata3D_mpi PRIVATE ${gfortran_opts})

add_library(neutraldata3D_geom_mpi OBJECT neutraldata3Dobj_geom_mpi.f90)
target_link_libraries(neutraldata3D_geom_mpi PRIVATE const config grid meshobj interp timeutils neutraldata neutraldata3D_mpi inputdata reader MPI::MPI_Fortran h5fortran::h5fortran)
target_compile_options(neutraldata3D_geom_mpi PRIVATE ${gfortran_opts})

add_library(neutraldata3D_geog_mpi OBJECT neutraldata3Dobj_geog_mpi.f90)
target_link_libraries(neutraldata3D_geog_mpi PRIVATE const config grid meshobj interp timeutils neutraldata neutraldata3D_mpi inputdata reader MPI::MPI_Fortran h5fortran::h5fortran)
target_compile_options(neutraldata3D_geog_mpi PRIVATE ${gfortran_opts})

