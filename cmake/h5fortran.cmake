if(hdf5) # needs to be run each configure for h5fortran_proj

include(${CMAKE_CURRENT_LIST_DIR}/hdf5.cmake)

include(ExternalProject)

ExternalProject_Add(h5fortran_proj
  GIT_REPOSITORY https://github.com/scivision/h5fortran.git
  GIT_TAG v1.3.1
  INSTALL_COMMAND ""  # disables the install step for the external project
)

ExternalProject_Get_Property(h5fortran_proj BINARY_DIR)
set(h5fortran_BINARY_DIR ${BINARY_DIR} CACHE PATH "path to h5fortran")

endif()