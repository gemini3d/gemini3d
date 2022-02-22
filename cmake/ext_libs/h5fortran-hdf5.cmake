if(hdf5)
  find_package(h5fortran CONFIG REQUIRED)
else(hdf5)
  add_library(h5fortran ${PROJECT_SOURCE_DIR}/src/vendor/h5fortran_dummy.f90)

  add_library(h5fortran::h5fortran INTERFACE IMPORTED)
  target_link_libraries(h5fortran::h5fortran INTERFACE h5fortran)

  install(TARGETS h5fortran EXPORT ${PROJECT_NAME}-targets)
endif(hdf5)
