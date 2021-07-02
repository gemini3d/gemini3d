set_directory_properties(PROPERTIES LABELS io)

add_library(logging OBJECT logging.f90)

# --- error handling -- dump data to disk
add_library(errors OBJECT error_stop.f90)
target_link_libraries(errors PRIVATE const h5fortran::h5fortran)

if(BUILD_TESTING)
  if(f2018assumed_rank)
    # easier to use, but needs latest compilers. For dev debugging.
    add_library(debug_dump OBJECT debug_dump.f90)
    target_link_libraries(debug_dump PRIVATE const h5fortran::h5fortran)
    set_target_properties(debug_dump PROPERTIES EXCLUDE_FROM_ALL true)
  endif(f2018assumed_rank)

  if(test_shaky)
    # this test sometimes causes a segfault on exit, which breaks CTest WILL_FAIL
    # this is a long-standing open issue with CMake et al,
    # so we by default omit this test.
    add_executable(test_errors test_errors.f90)
    target_link_libraries(test_errors PRIVATE errors const h5fortran::h5fortran)

    add_test(NAME gemini:stop:zero COMMAND $<TARGET_FILE:test_errors> root)
    set_tests_properties(gemini:stop:zero PROPERTIES
      WILL_FAIL true
      LABELS unit)
    # FIXME: add a follow-on test using Fixtures that read and verify expected file dump
  endif()

endif(BUILD_TESTING)

add_library(sanity_check OBJECT sanity_check.f90 $<TARGET_OBJECTS:errors>)
target_link_libraries(sanity_check PRIVATE const)

# --- pathlib
add_library(pathlib OBJECT pathlib.f90)

if(WIN32)
  target_sources(pathlib PRIVATE pathlib_windows.f90)
else()
  target_sources(pathlib PRIVATE pathlib_unix.f90)
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  target_sources(pathlib PRIVATE pathlib_intel.f90)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  target_sources(pathlib PRIVATE pathlib_gcc.f90)
else()
  message(WARNING "pathlib_dummy: better to write pathlib_gcc-like routine for compiler ${CMAKE_Fortran_COMPILER_ID}")
  target_sources(pathlib PRIVATE pathlib_dummy.f90)
endif()

if(BUILD_TESTING)
  add_executable(test_pathlib test_pathlib.f90)
  target_link_libraries(test_pathlib PRIVATE pathlib)
  add_test(NAME unit:pathlib COMMAND $<TARGET_FILE:test_pathlib>)
  set_tests_properties(unit:pathlib PROPERTIES
    LABELS unit
    TIMEOUT 10)
endif()

# --- config
add_library(config OBJECT config.f90 config_nml.f90 config_ini.f90)
target_link_libraries(config PRIVATE pathlib const)

set(_test_file ${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests/config/config_example.nml)

if(BUILD_TESTING)

add_executable(test_namelist test_nml.f90)
add_test(NAME unit:namelist COMMAND $<TARGET_FILE:test_namelist> ${_test_file})
set_tests_properties(unit:namelist PROPERTIES
REQUIRED_FILES ${_test_file}
LABELS unit
TIMEOUT 10)

add_executable(test_compiler_vendor test_compiler_vendor.f90)
target_link_libraries(test_compiler_vendor PRIVATE config pathlib)
add_test(NAME unit:compiler_vendor COMMAND $<TARGET_FILE:test_compiler_vendor>)
set_tests_properties(unit:compiler_vendor PROPERTIES
PASS_REGULAR_EXPRESSION ${CMAKE_Fortran_COMPILER_ID}
LABELS unit
TIMEOUT 10)

add_executable(test_expand_envvar test_expand_envvar.f90)
target_link_libraries(test_expand_envvar PRIVATE config pathlib)
add_test(NAME unit:expand_envvar COMMAND $<TARGET_FILE:test_expand_envvar>)
set_tests_properties(unit:expand_envvar PROPERTIES
ENVIRONMENT test__gem=hello
LABELS unit
TIMEOUT 10)

endif(BUILD_TESTING)

# --- io
configure_file(output.in.f90 output.f90 @ONLY)
add_library(io OBJECT io.f90 ${CMAKE_CURRENT_BINARY_DIR}/output.f90
  aurora.f90 aurora_raw.f90 aurora_hdf5.f90 aurora_ncdf.f90
  mag.f90 mag_hdf5.f90 mag_ncdf.f90 mag_raw.f90
  plasma_input.f90 plasma_input_hdf5.f90 plasma_input_ncdf.f90 plasma_input_raw.f90
  plasma_output.f90 plasma_output_hdf5.f90 plasma_output_ncdf.f90 plasma_output_raw.f90
  milestone.f90)

target_link_libraries(io PRIVATE calculus config const reader grid mpimod pathlib sanity_check timeutils
  h5fortran::h5fortran nc4fortran::nc4fortran)

# --- reader
add_library(reader OBJECT reader.f90 reader_hdf5.f90 reader_ncdf.f90 reader_raw.f90)
target_link_libraries(reader PRIVATE const pathlib h5fortran::h5fortran nc4fortran::nc4fortran)

if(netcdf AND BUILD_TESTING)
  add_executable(test_simsize_nc4 test_simsize_nc4.f90)
  target_link_libraries(test_simsize_nc4 PRIVATE nc4fortran::nc4fortran)
  # add_test(NAME unit:simsize_nc4 COMMAND test_simsize_nc4)
endif()

set_targ_props(config errors io logging pathlib reader sanity_check)
