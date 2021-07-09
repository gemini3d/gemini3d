set_directory_properties(PROPERTIES LABELS hwm)

add_library(hwm_ifc)
set_target_properties(hwm_ifc PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
target_include_directories(hwm_ifc INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>)

install(TARGETS hwm_ifc EXPORT ${PROJECT_NAME}Targets)

if(NOT hwm14)
  target_sources(hwm_ifc PRIVATE dummy.f90)
  return()
endif()

# --- build abstract HWM interface
# default is HWM14

include(${PROJECT_SOURCE_DIR}/cmake/ext_libs/hwm14.cmake)

target_sources(hwm_ifc PRIVATE hwm_interface.f90)
target_link_libraries(hwm_ifc PRIVATE HWM14::HWM14)

if(BUILD_TESTING)
  add_executable(test_hwm test_hwm.f90)
  target_link_libraries(test_hwm PRIVATE hwm_ifc)
  set_target_properties(test_hwm PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
    LABELS unit)

  add_test(NAME unit:hwm14 COMMAND $<TARGET_FILE:test_hwm>
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
  set_tests_properties(unit:hwm14 PROPERTIES
    TIMEOUT 5
    LABELS unit)
endif()
