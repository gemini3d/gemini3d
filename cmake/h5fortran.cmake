include(${CMAKE_CURRENT_LIST_DIR}/hdf5.cmake)

include(ExternalProject)

ExternalProject_Add(h5fortran_proj
  GIT_REPOSITORY https://github.com/scivision/h5fortran.git
  GIT_TAG v2.3.0
  INSTALL_COMMAND ""  # disables the install step for the external project
)

ExternalProject_Get_Property(h5fortran_proj BINARY_DIR)
#set(h5fortran_INCLUDE_DIR ${BINARY_DIR} CACHE PATH "path to h5fortran module files")

add_library(h5fortran STATIC IMPORTED GLOBAL)
add_dependencies(h5fortran h5fortran_proj)
set_target_properties(h5fortran PROPERTIES
  IMPORTED_LOCATION ${BINARY_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}h5fortran${CMAKE_STATIC_LIBRARY_SUFFIX}
  INTERFACE_INCLUDE_DIRECTORIES ${BINARY_DIR})
