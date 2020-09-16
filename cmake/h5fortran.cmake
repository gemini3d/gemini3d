if(hdf5)

  set(h5fortran_BUILD_TESTING false CACHE BOOL "h5fortran no test")

  if(MSVC AND NOT DEFINED HDF5OK)
  include(${CMAKE_CURRENT_LIST_DIR}/win32_hdf5.cmake)
  endif()

  find_package(h5fortran CONFIG)
  if(h5fortran_FOUND)
    include(${h5fortran_DIR}/h5fortranTargets.cmake)
  else()
    include(FetchContent)
    FetchContent_Declare(h5fortran_proj
      GIT_REPOSITORY https://github.com/geospace-code/h5fortran.git
      GIT_TAG v3.2.2)

    FetchContent_MakeAvailable(h5fortran_proj)
  endif()

  if(NOT HDF5OK)
    message(FATAL_ERROR "HDF5 was requested but is not available.")
  endif()

else(hdf5)
  message(VERBOSE " using h5fortran dummy")

  add_library(h5fortran ${CMAKE_CURRENT_SOURCE_DIR}/src/vendor/h5fortran/dummy.f90)
  target_include_directories(h5fortran INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
  set_target_properties(h5fortran PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
  add_library(h5fortran::h5fortran ALIAS h5fortran)
endif(hdf5)
