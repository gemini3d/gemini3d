# --- build abstract HWM interface
# default is HWM14

include(${PROJECT_SOURCE_DIR}/cmake/ext_libs/hwm14.cmake)

add_library(hwm_ifc OBJECT hwm_interface.f90)
target_link_libraries(hwm_ifc PRIVATE HWM14::HWM14)
target_include_directories(hwm_ifc INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(hwm_ifc PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

if(${PROJECT_NAME}_BUILD_TESTING)
  add_executable(test_hwm test_hwm.f90)
  target_link_libraries(test_hwm PRIVATE hwm_ifc)
  set_target_properties(test_hwm PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

  add_test(NAME unit:hwm14 COMMAND $<TARGET_FILE:test_hwm>
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
  set_tests_properties(unit:hwm14 PROPERTIES TIMEOUT 5)
endif()
