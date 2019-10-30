if(usehdf5)

include(ExternalProject)

ExternalProject_Add(h5fortran
  GIT_REPOSITORY https://github.com/scivision/h5fortran.git
  GIT_TAG master  # FIXME: it's better to use a specific Git revision or Git tag for reproducibility
  INSTALL_COMMAND ""  # disables the install step for the external project
)

ExternalProject_Get_Property(h5fortran BINARY_DIR)
set(h5fortran_BINARY_DIR BINARY_DIR)  # just to avoid accidentally reusing the variable name.

endif()