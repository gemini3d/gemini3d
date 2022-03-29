if(hdf5)
  find_package(h5fortran CONFIG REQUIRED)
else()
  include(${CMAKE_CURRENT_LIST_DIR}/../package/StubPackage.cmake)

  stub_package(h5fortran)

  add_library(h5fortran ${PROJECT_SOURCE_DIR}/src/vendor/h5fortran_dummy.f90)

  add_library(h5fortran::h5fortran INTERFACE IMPORTED)
  target_link_libraries(h5fortran::h5fortran INTERFACE h5fortran)

  install(TARGETS h5fortran EXPORT h5fortran-targets)
  install(FILES ${PROJECT_BINARY_DIR}/include/h5fortran.mod TYPE INCLUDE)
endif()
