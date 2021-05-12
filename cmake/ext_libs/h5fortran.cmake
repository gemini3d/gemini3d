include(ExternalProject)

if(hdf5)
  find_package(h5fortran CONFIG)
  if(h5fortran_FOUND)
    return()
  endif()

  if(NOT hdf5_external)
    find_package(HDF5 COMPONENTS Fortran HL)
  endif()

  if(NOT HDF5_FOUND OR hdf5_external)
    include(${CMAKE_CURRENT_LIST_DIR}/build_hdf5.cmake)
  endif()

  if(NOT h5fortran_ROOT)
    if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        set(h5fortran_ROOT ${PROJECT_BINARY_DIR})
      else()
        set(h5fortran_ROOT ${CMAKE_INSTALL_PREFIX})
    endif()
  endif()
  if(HDF5_FOUND)
    set(HDF5_ROOT ${HDF5_INCLUDE_DIR}/..)
  else()
    set(HDF5_ROOT ${h5fortran_ROOT})
  endif()

  message(VERBOSE "HDF5_ROOT: ${HDF5_ROOT}")

  set(h5fortran_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/include)
  set(h5fortran_LIBRARIES ${h5fortran_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}h5fortran${CMAKE_STATIC_LIBRARY_SUFFIX})

  ExternalProject_Add(H5FORTRAN
    GIT_REPOSITORY ${h5fortran_git}
    GIT_TAG ${h5fortran_tag}
    CMAKE_ARGS -DHDF5_ROOT:PATH=${HDF5_ROOT} -DCMAKE_INSTALL_PREFIX:PATH=${h5fortran_ROOT} -DBUILD_SHARED_LIBS:BOOL=false -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING:BOOL=false
    BUILD_BYPRODUCTS ${h5fortran_LIBRARIES}
    INACTIVITY_TIMEOUT 15
    CONFIGURE_HANDLED_BY_BUILD ON
    )

  file(MAKE_DIRECTORY ${h5fortran_INCLUDE_DIRS})

  add_library(h5fortran::h5fortran INTERFACE IMPORTED)
  target_link_libraries(h5fortran::h5fortran INTERFACE ${h5fortran_LIBRARIES})
  target_include_directories(h5fortran::h5fortran INTERFACE ${h5fortran_INCLUDE_DIRS})

  # race condition for linking without this
  add_dependencies(h5fortran::h5fortran H5FORTRAN)

  target_link_libraries(h5fortran::h5fortran INTERFACE HDF5::HDF5)

else(hdf5)
  message(VERBOSE " using h5fortran dummy")

  add_library(h5fortran ${CMAKE_CURRENT_SOURCE_DIR}/src/vendor/h5fortran_dummy.f90)
  target_include_directories(h5fortran INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
  set_target_properties(h5fortran PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
  add_library(h5fortran::h5fortran ALIAS h5fortran)
endif(hdf5)
