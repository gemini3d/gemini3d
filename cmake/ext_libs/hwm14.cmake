if(hwm14)
  find_package(hwm14 CONFIG REQUIRED)
else()
  include(${CMAKE_CURRENT_LIST_DIR}/../package/StubPackage.cmake)

  stub_package(hwm14)

  add_library(hwm_ifc ${CMAKE_CURRENT_SOURCE_DIR}/src/vendor/hwm14_dummy.f90)

  add_library(hwm14::hwm_ifc INTERFACE IMPORTED)
  target_link_libraries(hwm14::hwm_ifc INTERFACE hwm_ifc)

  install(TARGETS hwm_ifc EXPORT hwm14-targets)
  install(FILES ${PROJECT_BINARY_DIR}/include/hwm_interface.mod TYPE INCLUDE)
endif()
