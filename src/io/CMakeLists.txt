add_library(logging OBJECT logging.f90)

add_library(sanity_check OBJECT sanity_check.f90)
target_link_libraries(sanity_check PRIVATE const)

# --- pathlib
if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  set(dir_exist "inquire(directory=path, exist=exists)")
else()
  set(dir_exist "inquire(file=path, exist=exists)")
endif()
configure_file(pathlib.in.f90 pathlib.f90)

add_library(pathlib OBJECT ${CMAKE_CURRENT_BINARY_DIR}/pathlib.f90)
if(WIN32)
  target_sources(pathlib PRIVATE pathlib_windows.f90)
else()
  target_sources(pathlib PRIVATE pathlib_unix.f90)
endif()

add_executable(test_pathlib test_pathlib.f90)
target_link_libraries(test_pathlib pathlib)
add_test(NAME unit:pathlib COMMAND $<TARGET_FILE:test_pathlib>)

# --- config
add_library(config OBJECT config.f90 config_nml.f90 config_ini.f90)
target_link_libraries(config PRIVATE pathlib const)

set(_test_file ${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests/config/config_example.nml)

add_executable(test_namelist test_nml.f90)
add_test(NAME unit:namelist COMMAND $<TARGET_FILE:test_namelist> ${_test_file})
set_tests_properties(unit:namelist PROPERTIES
REQUIRED_FILES ${_test_file}
TIMEOUT 10)

add_executable(test_compiler_vendor test_compiler_vendor.f90)
target_link_libraries(test_compiler_vendor config pathlib)
add_test(NAME unit:compiler_vendor COMMAND $<TARGET_FILE:test_compiler_vendor>)
set_tests_properties(unit:compiler_vendor PROPERTIES PASS_REGULAR_EXPRESSION ${CMAKE_Fortran_COMPILER_ID})

# --- io
if(mpi)
add_library(io OBJECT io.f90 output.f90
  aurora.f90 aurora_raw.f90 aurora_hdf5.f90 aurora_ncdf.f90
  mag.f90 mag_hdf5.f90 mag_ncdf.f90 mag_raw.f90
  plasma_input.f90 plasma_input_hdf5.f90 plasma_input_ncdf.f90 plasma_input_raw.f90
  plasma_output.f90 plasma_output_hdf5.f90 plasma_output_ncdf.f90 plasma_output_raw.f90
  milestone.f90)

target_link_libraries(io PRIVATE calculus config const reader grid mpimod sanity_check timeutils
  h5fortran::h5fortran nc4fortran::nc4fortran)
endif(mpi)

# --- reader
add_library(reader OBJECT reader.f90 reader_hdf5.f90 reader_ncdf.f90 reader_raw.f90)
target_link_libraries(reader PRIVATE const pathlib h5fortran::h5fortran nc4fortran::nc4fortran)

if(netcdf)
  add_executable(test_simsize_nc4 test_simsize_nc4.f90)
  target_link_libraries(test_simsize_nc4 nc4fortran::nc4fortran)
  # add_test(NAME unit:simsize_nc4 COMMAND test_simsize_nc4)
endif(netcdf)
