if(netcdf)
include(FetchContent)

FetchContent_Declare(nc4fortran_proj
  GIT_REPOSITORY https://github.com/scivision/nc4fortran.git
  GIT_TAG v0.4.2
)

FetchContent_MakeAvailable(nc4fortran_proj)
endif()

if(NOT netcdf OR NOT NCDFOK)
  message(VERBOSE "NetCDF4 nc4fortran dummy library")

  add_library(nc4fortran ${CMAKE_CURRENT_SOURCE_DIR}/src/vendor/nc4fortran/dummy.f90)
  target_include_directories(nc4fortran INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
  set_target_properties(nc4fortran PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
  add_library(nc4fortran::nc4fortran ALIAS nc4fortran)
endif()