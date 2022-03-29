if(netcdf)
  find_package(nc4fortran CONFIG REQUIRED)
else()
  include(${CMAKE_CURRENT_LIST_DIR}/../package/StubPackage.cmake)

  stub_package(nc4fortran)

  add_library(nc4fortran ${CMAKE_CURRENT_SOURCE_DIR}/src/vendor/nc4fortran_dummy.f90)

  add_library(nc4fortran::nc4fortran INTERFACE IMPORTED)
  target_link_libraries(nc4fortran::nc4fortran INTERFACE nc4fortran)

  install(TARGETS nc4fortran EXPORT nc4fortran-targets)
  install(FILES ${PROJECT_BINARY_DIR}/include/nc4fortran.mod TYPE INCLUDE)
endif()
