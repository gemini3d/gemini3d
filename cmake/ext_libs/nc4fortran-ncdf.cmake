if(netcdf)
  include(${CMAKE_CURRENT_LIST_DIR}/netcdf.cmake)
else(netcdf)
  add_library(nc4fortran ${CMAKE_CURRENT_SOURCE_DIR}/src/vendor/nc4fortran_dummy.f90)

  add_library(nc4fortran::nc4fortran INTERFACE IMPORTED)
  target_link_libraries(nc4fortran::nc4fortran INTERFACE nc4fortran)

  install(TARGETS nc4fortran EXPORT ${PROJECT_NAME}-targets)
endif(netcdf)
