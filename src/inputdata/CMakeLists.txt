set_directory_properties(PROPERTIES LABELS inputdata)

add_library(inputdata OBJECT inputdataobj.f90)
target_link_libraries(inputdata PRIVATE const config meshobj interp timeutils)

add_library(precipdata OBJECT precipdataobj.f90)
target_link_libraries(precipdata PRIVATE const config meshobj interp timeutils mpimod inputdata MPI::MPI_Fortran)

add_library(efielddata OBJECT efielddataobj.f90)
target_link_libraries(efielddata PRIVATE const config grid meshobj interp timeutils inputdata)

add_library(neutraldata OBJECT neutraldataobj.f90)
target_link_libraries(neutraldata PRIVATE const inputdata)

add_library(neutraldata3D OBJECT neutraldata3Dobj.f90 neuslab.f90)
target_link_libraries(neutraldata3D PRIVATE const config grid meshobj interp timeutils neutraldata inputdata MPI::MPI_Fortran h5fortran::h5fortran)
target_compile_options(neutraldata3D PRIVATE ${gfortran_opts})
