if(usehdf5)

include(${PROJECT_SOURCE_DIR}/cmake/hdf5.cmake)

include(ExternalProject)

ExternalProject_Add(h5fortran_proj
  GIT_REPOSITORY https://github.com/scivision/h5fortran.git
  GIT_TAG v1.3.0
  INSTALL_COMMAND ""  # disables the install step for the external project
)

ExternalProject_Get_Property(h5fortran_proj BINARY_DIR)
set(h5fortran_BINARY_DIR ${BINARY_DIR})  # just to avoid accidentally reusing the variable name.

message(STATUS "Using h5fortran in ${h5fortran_BINARY_DIR}")

endif()