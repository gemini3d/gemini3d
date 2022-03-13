if(netcdf)
  find_package(nc4fortran CONFIG REQUIRED)
else(netcdf)
  # also need package for this dummy
  include(CMakePackageConfigHelpers)

  configure_package_config_file(${PROJECT_SOURCE_DIR}/cmake/package/nc4fortran-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/nc4fortran-config.cmake
  INSTALL_DESTINATION cmake
  )

  write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/nc4fortran-config-version.cmake
  COMPATIBILITY SameMinorVersion
  )

  install(EXPORT nc4fortran-targets
  NAMESPACE nc4fortran::
  DESTINATION cmake
  )

  install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/nc4fortran-config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/nc4fortran-config-version.cmake
  DESTINATION cmake
  )

  add_library(nc4fortran ${CMAKE_CURRENT_SOURCE_DIR}/src/vendor/nc4fortran_dummy.f90)

  add_library(nc4fortran::nc4fortran INTERFACE IMPORTED)
  target_link_libraries(nc4fortran::nc4fortran INTERFACE nc4fortran)

  install(TARGETS nc4fortran EXPORT nc4fortran-targets)
  install(FILES ${PROJECT_BINARY_DIR}/include/nc4fortran.mod TYPE INCLUDE)
endif(netcdf)
