include(${CMAKE_CURRENT_LIST_DIR}/netcdf.cmake)

include(ExternalProject)

ExternalProject_Add(nc4fortran_proj
  GIT_REPOSITORY https://github.com/scivision/nc4fortran.git
  GIT_TAG master
  INSTALL_COMMAND ""  # disables the install step for the external project
)

ExternalProject_Get_Property(nc4fortran_proj BINARY_DIR)

add_library(nc4fortran STATIC IMPORTED GLOBAL)
add_dependencies(nc4fortran nc4fortran_proj)
set_target_properties(nc4fortran PROPERTIES
  IMPORTED_LOCATION ${BINARY_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}nc4fortran${CMAKE_STATIC_LIBRARY_SUFFIX}
  INTERFACE_INCLUDE_DIRECTORIES ${BINARY_DIR})
