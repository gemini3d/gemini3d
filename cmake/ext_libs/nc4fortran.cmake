# Leave nc4fortran as FetchContent as we use wrangle NetCDF library distinctions there
include(ExternalProject)

if(netcdf)
  find_package(nc4fortran CONFIG)
  if(nc4fortran_FOUND)
    return()
  endif()

  find_package(NetCDF REQUIRED COMPONENTS Fortran)

  if(NOT nc4fortran_ROOT)
    if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        set(nc4fortran_ROOT ${PROJECT_BINARY_DIR})
      else()
        set(nc4fortran_ROOT ${CMAKE_INSTALL_PREFIX})
    endif()
  endif()

  if(NOT DEFINED NetCDF_ROOT)
    set(NetCDF_ROOT ${NetCDF_C_INCLUDE_DIRS}/..)
  endif()
  message(VERBOSE "NetCDF_ROOT: ${NetCDF_ROOT}")

  set(nc4fortran_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/include)
  set(nc4fortran_LIBRARIES ${nc4fortran_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}nc4fortran${CMAKE_STATIC_LIBRARY_SUFFIX})

  ExternalProject_Add(NC4FORTRAN
    GIT_REPOSITORY ${nc4fortran_git}
    GIT_TAG ${nc4fortran_tag}
    CMAKE_ARGS -DNetCDF_ROOT:PATH=${NetCDF_ROOT} -DCMAKE_INSTALL_PREFIX:PATH=${nc4fortran_ROOT} -DBUILD_SHARED_LIBS:BOOL=false -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING:BOOL=false
    BUILD_BYPRODUCTS ${nc4fortran_LIBRARIES}
    INACTIVITY_TIMEOUT 30
    CONFIGURE_HANDLED_BY_BUILD ON
    )

  file(MAKE_DIRECTORY ${nc4fortran_INCLUDE_DIRS})

  add_library(nc4fortran::nc4fortran INTERFACE IMPORTED)
  target_link_libraries(nc4fortran::nc4fortran INTERFACE ${nc4fortran_LIBRARIES})
  target_include_directories(nc4fortran::nc4fortran INTERFACE ${nc4fortran_INCLUDE_DIRS})

  # race condition for linking without this
  add_dependencies(nc4fortran::nc4fortran NC4FORTRAN)

  target_link_libraries(nc4fortran::nc4fortran INTERFACE NetCDF::NetCDF_Fortran)

else(netcdf)
  message(VERBOSE " using nc4fortran dummy")

  add_library(nc4fortran ${CMAKE_CURRENT_SOURCE_DIR}/src/vendor/nc4fortran_dummy.f90)
  target_include_directories(nc4fortran INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
  set_target_properties(nc4fortran PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
  add_library(nc4fortran::nc4fortran ALIAS nc4fortran)
endif(netcdf)
