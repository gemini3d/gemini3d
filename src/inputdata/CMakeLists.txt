set_directory_properties(PROPERTIES LABELS inputdata)

add_library(inputdata OBJECT inputdataobj.f90)
target_link_libraries(inputdata PRIVATE const config meshobj interp timeutils)

add_library(precipdata OBJECT precipdataobj.f90)
target_link_libraries(precipdata PRIVATE const config meshobj interp timeutils mpimod inputdata MPI::MPI_Fortran)

add_library(efielddata OBJECT efielddataobj.f90)
target_link_libraries(efielddata PRIVE const config meshobj interp timeutils inputdata)
