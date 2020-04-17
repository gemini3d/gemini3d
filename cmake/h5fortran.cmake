if(hdf5)
  include(FetchContent)

  FetchContent_Declare(h5fortran_proj
    GIT_REPOSITORY https://github.com/scivision/h5fortran.git
    GIT_TAG v2.9.1)

  FetchContent_MakeAvailable(h5fortran_proj)

  if(NOT HDF5OK)
    message(FATAL_ERROR "HDF5 was requested but is not available.")
  endif()
else()
  message(VERBOSE " using HDF5 h5fortran dummy library")

  add_library(h5fortran ${CMAKE_CURRENT_SOURCE_DIR}/src/vendor/h5fortran/dummy.f90)
  target_include_directories(h5fortran INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
  set_target_properties(h5fortran PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
  add_library(h5fortran::h5fortran ALIAS h5fortran)
endif()