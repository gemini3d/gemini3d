if(netcdf)

  if(CMAKE_VERSION VERSION_LESS 3.19)
    message(FATAL_ERROR "NetCDF4 requires CMake >= 3.19")
  endif()

  set(nc4fortran_BUILD_TESTING false CACHE BOOL "h5fortran no test")

  find_package(nc4fortran CONFIG)
  if(nc4fortran_FOUND)
    include(${nc4fortran_DIR}/nc4fortranTargets.cmake)
  else()
    include(FetchContent)
    FetchContent_Declare(ncproj
      GIT_REPOSITORY ${nc4fortran_git}
      GIT_TAG ${nc4fortran_tag}
      GIT_SHALLOW true)

    FetchContent_MakeAvailable(ncproj)
  endif()

  if(NOT TARGET NetCDF::NetCDF_Fortran)
    message(FATAL_ERROR "NetCDF4 was requested but is not available.")
  endif()

  set(NetCDF_FOUND true)
else(netcdf)
  message(VERBOSE " using nc4fortran dummy")

  add_library(nc4fortran ${CMAKE_CURRENT_SOURCE_DIR}/src/vendor/nc4fortran/dummy.f90)
  target_include_directories(nc4fortran INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
  set_target_properties(nc4fortran PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
  add_library(nc4fortran::nc4fortran ALIAS nc4fortran)

  set(NetCDF_FOUND false)
endif(netcdf)
