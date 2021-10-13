include(ExternalProject)

if(hdf5)

  if(NOT hdf5_external)
    # h5fortran inside if() because h5fortran config calls find_package(HDF5)
    # disabled h5fortran search for now because we're undergoing rapid devel soon for MPI.
    # find_package(h5fortran CONFIG)
    # if(h5fortran_FOUND)
    #   return()
    # endif()

    if(autobuild)
      find_package(HDF5 COMPONENTS Fortran)
    else()
      find_package(HDF5 COMPONENTS Fortran REQUIRED)
    endif()
  endif()

  if(NOT HDF5_FOUND OR hdf5_external)
    include(${CMAKE_CURRENT_LIST_DIR}/build_hdf5.cmake)
  else()
    add_custom_target(HDF5)
  endif()

  if(NOT h5fortran_ROOT)
    if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
      set(h5fortran_ROOT ${PROJECT_BINARY_DIR} CACHE PATH "default ROOT")
    else()
      set(h5fortran_ROOT ${CMAKE_INSTALL_PREFIX})
    endif()
  endif()

  find_package(ZLIB)

  set(h5fortran_INCLUDE_DIRS ${h5fortran_ROOT}/include)
  set(h5fortran_LIBRARIES ${h5fortran_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}h5fortran${CMAKE_STATIC_LIBRARY_SUFFIX})

  set(h5fortran_cmake_args
  -DZLIB_ROOT:PATH=${ZLIB_ROOT}
  -DCMAKE_INSTALL_PREFIX:PATH=${h5fortran_ROOT}
  -DBUILD_SHARED_LIBS:BOOL=false
  -DCMAKE_BUILD_TYPE=Release
  -DBUILD_TESTING:BOOL=false
  -Dautobuild:BOOL=false
  )

  if(HDF5_ROOT)
    list(APPEND h5fortran_cmake_args -DHDF5_ROOT:PATH=${HDF5_ROOT})
  endif()

  ExternalProject_Add(H5FORTRAN
    GIT_REPOSITORY ${h5fortran_git}
    GIT_TAG ${h5fortran_tag}
    CMAKE_ARGS ${h5fortran_cmake_args}
    CMAKE_GENERATOR ${EXTPROJ_GENERATOR}
    BUILD_BYPRODUCTS ${h5fortran_LIBRARIES}
    INACTIVITY_TIMEOUT 15
    CONFIGURE_HANDLED_BY_BUILD ON
    DEPENDS HDF5)

  file(MAKE_DIRECTORY ${h5fortran_INCLUDE_DIRS})

  add_library(h5fortran::h5fortran INTERFACE IMPORTED)
  target_link_libraries(h5fortran::h5fortran INTERFACE ${h5fortran_LIBRARIES})
  target_include_directories(h5fortran::h5fortran INTERFACE ${h5fortran_INCLUDE_DIRS})

  # race condition for linking without this
  add_dependencies(h5fortran::h5fortran H5FORTRAN)

  target_link_libraries(h5fortran::h5fortran INTERFACE HDF5::HDF5)

else(hdf5)
  message(VERBOSE "using h5fortran dummy")

  add_library(h5fortran ${CMAKE_CURRENT_SOURCE_DIR}/src/vendor/h5fortran_dummy.f90)
  add_library(h5fortran::h5fortran ALIAS h5fortran)
endif(hdf5)
