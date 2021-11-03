set_directory_properties(PROPERTIES LABELS msis)

# --- MSISE00 and MSIS 2.0 require legacy workaround due to non-standard Fortran code

# Gfortran >= 8 need -Wno-pedantic to allow mismatched array size inhernet to MSIS.
# "-w" doesn't disable pedantic
set(msis00_gnu_flags -w -std=legacy -Wno-pedantic -fno-implicit-none -Wno-error=array-bounds -fcheck=no-all)
# msis2.0 will intermittantly fail with -std=legacy with all zero outputs, e.g. GCC 8.5.0
set(msis2_gnu_flags -w -Wno-error=array-bounds -fcheck=no-all)

set(msis_intel_flags -nowarn)

# --- MSISE00 is always available, but MSIS 2.0 might not be available.

add_library(msis00mod msis00_gfortran.f msise00_data.f)
# keep msis00 NOT an OBJECT to avoid having to link it explicitly everywhere
target_compile_options(msis00mod PRIVATE
"$<$<Fortran_COMPILER_ID:GNU>:${msis00_gnu_flags}>"
"$<$<Fortran_COMPILER_ID:Intel,IntelLLVM>:${msis_intel_flags}>")

# --- build abstract MSIS interface

add_library(msis_ifc OBJECT msis_interface.f90)
target_link_libraries(msis_ifc PRIVATE msis00mod)

# --- optional MSIS 2.0

if(msis2)
  include(${PROJECT_SOURCE_DIR}/cmake/ext_libs/msis2.cmake)
  target_link_libraries(msis_ifc PRIVATE msis2)
  target_compile_options(msis2 PRIVATE
  "$<$<Fortran_COMPILER_ID:GNU>:${msis2_gnu_flags}>"
  "$<$<Fortran_COMPILER_ID:Intel,IntelLLVM>:${msis_intel_flags}>")

  install(TARGETS msis2
    EXPORT ${PROJECT_NAME}Targets
  )
else(msis2)
  target_sources(msis_ifc PRIVATE msis2_dummy.f90)
endif(msis2)

# cmake bug? Object library not a good idea for msis 2.0 as have to specify everywhere
# install(TARGETS msis_ifc msis00mod EXPORT ${PROJECT_NAME}Targets)

if(NOT hdf5)
  return()
endif()

# --- for setting up an equilibrium simulation

add_executable(msis_setup msis_driver.f90)
target_link_libraries(msis_setup PRIVATE msis_ifc const h5fortran::h5fortran
$<$<BOOL:${HDF5_HAVE_PARALLEL}>:MPI::MPI_Fortran>
)
set_target_properties(msis_setup PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

install(TARGETS msis00mod msis_setup EXPORT ${PROJECT_NAME}Targets)

if(BUILD_TESTING)

add_executable(msis_generate test_gen_in.f90)
target_link_libraries(msis_generate PRIVATE h5fortran::h5fortran
$<$<BOOL:${HDF5_HAVE_PARALLEL}>:MPI::MPI_Fortran>
)

add_executable(msis_compare test_msis.f90)
target_link_libraries(msis_compare PRIVATE msis_ifc assert h5fortran::h5fortran
$<$<BOOL:${HDF5_HAVE_PARALLEL}>:MPI::MPI_Fortran>
)

# --- MSISE00 tests
set(msis00_in ${CMAKE_CURRENT_BINARY_DIR}/msis00_test_in.h5)
set(msis00_out ${CMAKE_CURRENT_BINARY_DIR}/msis00_test_out.h5)
set(msis00_ref ${PROJECT_SOURCE_DIR}/test_data/msis/msis00_test_out.h5)

add_test(NAME unit:msis00_generate
  COMMAND $<TARGET_FILE:msis_generate> ${msis00_in} 0)
set_tests_properties(unit:msis00_generate PROPERTIES
  FIXTURES_SETUP msis00_gen_fxt
  LABELS unit
  TIMEOUT 10
)

add_test(NAME unit:msis00_setup
  COMMAND $<TARGET_FILE:msis_setup> ${msis00_in} ${msis00_out})
set_tests_properties(unit:msis00_setup PROPERTIES
  FIXTURES_SETUP msis00_fxt
  FIXTURES_REQUIRED msis00_gen_fxt
  LABELS unit
  REQUIRED_FILES ${msis00_in}
  TIMEOUT 10
)

add_test(NAME unit:msis00_compare
  COMMAND $<TARGET_FILE:msis_compare> ${msis00_out} ${msis00_ref})
set_tests_properties(unit:msis00_compare PROPERTIES
  FIXTURES_REQUIRED msis00_fxt
  LABELS unit
  TIMEOUT 10
  REQUIRED_FILES "${msis00_out};${msis00_ref}"
)

set_tests_properties(unit:msis00_setup unit:msis00_compare PROPERTIES
ENVIRONMENT $<$<BOOL:${test_dll_path}>:"PATH=${test_dll_path}">
)

# --- MSIS 2.0 tests
if(msis2)

  set(msis2_in ${CMAKE_CURRENT_BINARY_DIR}/msis2_test_in.h5)
  set(msis2_out ${CMAKE_CURRENT_BINARY_DIR}/msis2_test_out.h5)
  set(msis2_ref ${PROJECT_SOURCE_DIR}/test_data/msis/msis2_test_out.h5)

  add_test(NAME unit:msis2_generate
    COMMAND $<TARGET_FILE:msis_generate> ${msis2_in} 20
  )
  set_tests_properties(unit:msis2_generate PROPERTIES
    FIXTURES_SETUP msis2_gen_fxt
    LABELS unit
    TIMEOUT 10
  )

  add_test(NAME unit:msis2_setup
    COMMAND $<TARGET_FILE:msis_setup> ${msis2_in} ${msis2_out}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )
  set_tests_properties(unit:msis2_setup PROPERTIES
    FIXTURES_SETUP msis2_fxt
    FIXTURES_REQUIRED msis2_gen_fxt
    LABELS unit
    REQUIRED_FILES ${msis2_in}
    TIMEOUT 10
  )

  add_test(NAME unit:msis2_compare
    COMMAND $<TARGET_FILE:msis_compare> ${msis2_out} ${msis2_ref}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )
  set_tests_properties(unit:msis2_compare PROPERTIES
    FIXTURES_REQUIRED msis2_fxt
    LABELS unit
    TIMEOUT 10
    REQUIRED_FILES "${msis2_out};${msis2_ref}"
  )

 set_tests_properties(unit:msis2_setup unit:msis2_compare PROPERTIES
 ENVIRONMENT $<$<BOOL:${test_dll_path}>:"PATH=${test_dll_path}">
 )

endif(msis2)

endif(BUILD_TESTING)
