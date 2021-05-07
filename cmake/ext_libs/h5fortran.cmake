include(ExternalProject)

if(hdf5)
  find_package(h5fortran CONFIG)
  if(h5fortran_FOUND)
    return()
  endif()

  if(NOT h5fortran_ROOT)
    if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        set(h5fortran_ROOT ${PROJECT_BINARY_DIR})
      else()
        set(h5fortran_ROOT ${CMAKE_INSTALL_PREFIX})
    endif()
  endif()

  set(h5fortran_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/include)
  set(h5fortran_LIBRARIES ${h5fortran_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}h5fortran${CMAKE_STATIC_LIBRARY_SUFFIX})

  ExternalProject_Add(H5FORTRAN
    GIT_REPOSITORY ${h5fortran_git}
    GIT_TAG ${h5fortran_tag}
    UPDATE_DISCONNECTED ${EP_UPDATE_DISCONNECTED}
    CMAKE_ARGS -Dhdf5_external:BOOL=${hdf5_external} -DHDF5_ROOT:PATH=${HDF5_ROOT} -DCMAKE_INSTALL_PREFIX:PATH=${h5fortran_ROOT} -DBUILD_SHARED_LIBS:BOOL=false -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING:BOOL=false
    BUILD_BYPRODUCTS ${h5fortran_LIBRARIES}
    INACTIVITY_TIMEOUT 30
    CONFIGURE_HANDLED_BY_BUILD ON
    )

  file(MAKE_DIRECTORY ${h5fortran_INCLUDE_DIRS})

  add_library(h5fortran::h5fortran INTERFACE IMPORTED)
  target_link_libraries(h5fortran::h5fortran INTERFACE ${h5fortran_LIBRARIES})
  target_include_directories(h5fortran::h5fortran INTERFACE ${h5fortran_INCLUDE_DIRS})

  # race condition for linking without this
  add_dependencies(h5fortran::h5fortran H5FORTRAN)

  # --- superproject workaround: HDF5
  if(NOT hdf5_external)
    find_package(HDF5 COMPONENTS Fortran HL)
  endif()
  # ExternalProject does not carry the CMake package configuration along.
  # so if h5fortran manually or auto builds HDF5, we need to point to that manually.
  # The canonical way to do this would be via a CMake superproject, but let's see if we
  # can avoid that CMake code duplication for now.
  if(NOT HDF5_FOUND)
    # we assume that h5fortran determined it needed to build HDF5.
    # a more canonical way would be to do like we do Mumps, Lapack and Scalapack--
    # a superproject like approach instead would be more robust.

    # the code below heavily follows h5fortran/cmake/build_{zlib,hdf5}.cmake
    set(HDF5_LIBRARIES)
    foreach(_name hdf5_hl_fortran hdf5_hl_f90cstub hdf5_fortran hdf5_f90cstub hdf5_hl hdf5)
      list(APPEND HDF5_LIBRARIES ${h5fortran_ROOT}/lib/lib${_name}${CMAKE_STATIC_LIBRARY_SUFFIX})
    endforeach()

    set(HDF5_INCLUDE_DIRS ${h5fortran_ROOT}/include)

    add_library(HDF5::HDF5 INTERFACE IMPORTED GLOBAL)
    target_include_directories(HDF5::HDF5 INTERFACE "${HDF5_INCLUDE_DIRS}")
    target_link_libraries(HDF5::HDF5 INTERFACE "${HDF5_LIBRARIES}")

    if(NOT TARGET ZLIB::ZLIB)
      if(WIN32)
        set(ZLIB_name ${CMAKE_STATIC_LIBRARY_PREFIX}zlibstatic${CMAKE_STATIC_LIBRARY_SUFFIX})
      else()
        set(ZLIB_name ${CMAKE_STATIC_LIBRARY_PREFIX}z${CMAKE_STATIC_LIBRARY_SUFFIX})
      endif()

      set(ZLIB_INCLUDE_DIR ${h5fortran_ROOT}/include)
      set(ZLIB_LIBRARY ${h5fortran_ROOT}/lib/${ZLIB_name})

      add_library(ZLIB::ZLIB INTERFACE IMPORTED)
      add_dependencies(ZLIB::ZLIB ZLIB)  # to avoid include directory race condition
      target_link_libraries(ZLIB::ZLIB INTERFACE ${ZLIB_LIBRARY})
      target_include_directories(ZLIB::ZLIB INTERFACE ${ZLIB_INCLUDE_DIR})
    endif()

    target_link_libraries(HDF5::HDF5 INTERFACE ZLIB::ZLIB)

    set(THREADS_PREFER_PTHREAD_FLAG true)
    find_package(Threads)
    if(Threads_FOUND)
      target_link_libraries(HDF5::HDF5 INTERFACE Threads::Threads)
    endif(Threads_FOUND)

    # libdl and libm are needed on some systems--don't remove
    target_link_libraries(HDF5::HDF5 INTERFACE ${CMAKE_DL_LIBS})

    if(UNIX)
      target_link_libraries(HDF5::HDF5 INTERFACE m)
    endif(UNIX)

  endif()

  target_link_libraries(h5fortran::h5fortran INTERFACE HDF5::HDF5)
  # --- end superproject workaround

else(hdf5)
  message(VERBOSE " using h5fortran dummy")

  add_library(h5fortran ${CMAKE_CURRENT_SOURCE_DIR}/src/vendor/h5fortran_dummy.f90)
  target_include_directories(h5fortran INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
  set_target_properties(h5fortran PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
  add_library(h5fortran::h5fortran ALIAS h5fortran)
endif(hdf5)
