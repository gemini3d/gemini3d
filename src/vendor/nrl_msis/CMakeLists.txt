set_directory_properties(PROPERTIES LABELS msis)

# We assume MSISE00 is always available, but MSIS 2.0 might not be available.

add_library(msis00mod msis00_gfortran.f msise00_data.f)
# keep msis00 NOT an OBJECT to avoid having to link it explicitly everywhere
target_compile_options(msis00mod PRIVATE ${msis_flags})

# --- build abstract MSIS interface

add_library(msis_ifc OBJECT msis_interface.f90)
target_link_libraries(msis_ifc PRIVATE msis00mod)

# --- optional MSIS 2.0

if(msis20)
  include(${PROJECT_SOURCE_DIR}/cmake/ext_libs/msis2.cmake)
  target_link_libraries(msis_ifc PRIVATE msis2)
  target_compile_options(msis2 PRIVATE ${msis_flags})
else(msis20)
  target_sources(msis_ifc PRIVATE msis2_dummy.f90)
endif(msis20)

# cmake bug? Object library not a good idea for msis 2.0 as have to specify everywhere
# install(TARGETS msis_ifc msis00mod EXPORT ${PROJECT_NAME}Targets)

if(NOT hdf5)
  return()
endif()

# --- for setting up an equilibrium simulation

add_executable(msis_setup msis_driver.f90)
target_link_libraries(msis_setup PRIVATE msis_ifc const h5fortran::h5fortran)
set_target_properties(msis_setup PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
if(HDF5_HAVE_PARALLEL)
  target_link_libraries(msis_setup PRIVATE MPI::MPI_Fortran)
endif()

install(TARGETS msis_setup EXPORT ${PROJECT_NAME}Targets)

if(BUILD_TESTING)

set(msis_in ${CMAKE_CURRENT_BINARY_DIR}/msis_test_in.h5)
set(msis_out ${CMAKE_CURRENT_BINARY_DIR}/msis_test_out.h5)
set(msis_ref ${PROJECT_SOURCE_DIR}/test_data/msis_test_out.h5)

add_executable(msis_generate test_gen_in.f90)
target_link_libraries(msis_generate PRIVATE h5fortran::h5fortran)
if(HDF5_HAVE_PARALLEL)
  target_link_libraries(msis_generate PRIVATE MPI::MPI_Fortran)
endif()

add_test(NAME unit:msis_generate
  COMMAND $<TARGET_FILE:msis_generate> ${msis_in})
set_tests_properties(unit:msis_generate PROPERTIES
  FIXTURES_SETUP msis_gen_fxt
  LABELS unit
  TIMEOUT 10)

add_test(NAME unit:msis_setup
  COMMAND $<TARGET_FILE:msis_setup> ${msis_in} ${msis_out} 0)
set_tests_properties(unit:msis_setup PROPERTIES
  FIXTURES_SETUP msis_fxt
  FIXTURES_REQUIRED msis_gen_fxt
  LABELS unit
  REQUIRED_FILES ${msis_in}
  TIMEOUT 10)

add_executable(msis_compare test_msis.f90)
target_link_libraries(msis_compare PRIVATE msis_ifc assert h5fortran::h5fortran)
if(HDF5_HAVE_PARALLEL)
  target_link_libraries(msis_compare PRIVATE MPI::MPI_Fortran)
endif()

add_test(NAME unit:msis_compare
  COMMAND $<TARGET_FILE:msis_compare> ${msis_out} ${msis_ref})
set_tests_properties(unit:msis_compare PROPERTIES
  FIXTURES_REQUIRED msis_fxt
  LABELS unit
  TIMEOUT 10
  REQUIRED_FILES "${msis_out};${msis_ref}")

if(test_dll_path)
  set_tests_properties(unit:msis_setup unit:msis_compare PROPERTIES
    ENVIRONMENT "PATH=${test_dll_path}")
endif()

endif()
