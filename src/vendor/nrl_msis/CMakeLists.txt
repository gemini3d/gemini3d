# We assume MSISE00 is always available, but MSIS 2.0 might not be available.

add_library(msis00mod msis00_gfortran.f msise00_data.f)
# keep msis00 NOT an OBJECT to avoid having to link it explicitly everywhere
target_compile_options(msis00mod PRIVATE ${msis_flags})
target_include_directories(msis00mod INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(msis00mod PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)


# --- build abstract MSIS interface

add_library(msis_ifc msis_interface.f90)
# msis_ifc not OBJECT to avoid linking issues on CMake 3.13
target_link_libraries(msis_ifc PRIVATE msis00mod)
target_include_directories(msis_ifc INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(msis_ifc PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

# --- optional MSIS 2.0

if(msis20)
  include(${PROJECT_SOURCE_DIR}/cmake/msis2.cmake)
  target_link_libraries(msis_ifc PRIVATE msis2)
  target_compile_options(msis2 PRIVATE ${msis_flags})
else(msis20)
  target_sources(msis_ifc PRIVATE msis2_dummy.f90)
endif(msis20)

# --- for setting up an equilibrium simulation

add_executable(msis_setup msis_driver.f90)
target_link_libraries(msis_setup PRIVATE msis_ifc const h5fortran::h5fortran)
target_compile_options(msis_setup PRIVATE ${static_flag})
set_target_properties(msis_setup PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
