set_directory_properties(PROPERTIES LABELS msis)

# We assume MSISE00 is always available, but MSIS 2.0 might not be available.

add_library(msis00mod msis00_gfortran.f msise00_data.f)
# keep msis00 NOT an OBJECT to avoid having to link it explicitly everywhere
target_compile_options(msis00mod PRIVATE ${msis_flags})
target_include_directories(msis00mod INTERFACE
$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
$<INSTALL_INTERFACE:include>)
set_target_properties(msis00mod PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)


# --- build abstract MSIS interface

add_library(msis_ifc msis_interface.f90)
# msis_ifc not OBJECT to avoid linking issues on CMake 3.13
target_link_libraries(msis_ifc PRIVATE msis00mod)
target_include_directories(msis_ifc INTERFACE
$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
$<INSTALL_INTERFACE:include>)
set_target_properties(msis_ifc PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

# --- optional MSIS 2.0

if(msis20)
  include(${PROJECT_SOURCE_DIR}/cmake/ext_libs/msis2.cmake)
  target_link_libraries(msis_ifc PRIVATE msis2)
  target_compile_options(msis2 PRIVATE ${msis_flags})
else(msis20)
  target_sources(msis_ifc PRIVATE msis2_dummy.f90)
endif(msis20)

install(TARGETS msis_ifc msis00mod EXPORT ${PROJECT_NAME}Targets)

# --- for setting up an equilibrium simulation

add_executable(msis_setup msis_driver.f90)
target_link_libraries(msis_setup PRIVATE msis_ifc const h5fortran::h5fortran)
set_target_properties(msis_setup PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

install(TARGETS msis_setup
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME DESTINATION bin)

if(BUILD_TESTING)

set(_in ${PROJECT_SOURCE_DIR}/test_data/msis_test_in.h5)
set(_out ${CMAKE_CURRENT_BINARY_DIR}/msis_test_out.h5)
set(_ref ${PROJECT_SOURCE_DIR}/test_data/msis_test_out.h5)

add_test(NAME unit:msis_setup
  COMMAND $<TARGET_FILE:msis_setup> ${_in} ${_out} 0)
set_tests_properties(unit:msis_setup PROPERTIES
  FIXTURES_SETUP msis_fxt
  LABELS unit
  REQUIRED_FILES ${_in}
  TIMEOUT 10)

add_executable(msis_compare test_msis.f90)
target_link_libraries(msis_compare PRIVATE msis_ifc assert h5fortran::h5fortran)

add_test(NAME unit:msis_compare
  COMMAND $<TARGET_FILE:msis_compare> ${_out} ${_ref})
set_tests_properties(unit:msis_compare PROPERTIES
  FIXTURES_REQUIRED msis_fxt
  LABELS unit
  TIMEOUT 10
  REQUIRED_FILES "${_out};${_ref}")

if(test_dll_path)
  set_tests_properties(unit:msis_setup unit:msis_compare PROPERTIES
    ENVIRONMENT "PATH=${test_dll_path}")
endif()

endif()
