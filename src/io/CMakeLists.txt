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


# --- config
add_library(config OBJECT config.f90 config_nml.f90 config_ini.f90)
target_link_libraries(config PRIVATE pathlib const)

set(_test_file ${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests/config/test2dew_glow/config.nml)

add_executable(test_namelist test_nml.f90)
add_test(NAME unit:namelist COMMAND $<TARGET_FILE:test_namelist> ${_test_file})
set_tests_properties(unit:namelist PROPERTIES
REQUIRED_FILES ${_test_file}
TIMEOUT 10)

add_executable(test_config test_config.f90)
target_link_libraries(test_config PRIVATE config const pathlib)
add_test(NAME unit:config COMMAND $<TARGET_FILE:test_config> ${_test_file})
set_tests_properties(unit:config PROPERTIES
REQUIRED_FILES ${_test_file}
TIMEOUT 10)

# --- io
if(MPI_OK)
  add_library(io OBJECT io.f90 output.f90
    aurora.f90 aurora_raw.f90 aurora_hdf5.f90 aurora_ncdf.f90
    mag.f90 mag_hdf5.f90 mag_ncdf.f90 mag_raw.f90
    plasma_input.f90 plasma_input_hdf5.f90 plasma_input_ncdf.f90 plasma_input_raw.f90
    plasma_output.f90 plasma_output_hdf5.f90 plasma_output_ncdf.f90 plasma_output_raw.f90)

  target_link_libraries(io PRIVATE calculus config const reader grid mpimod sanity_check timeutils
    h5fortran::h5fortran nc4fortran::nc4fortran)
endif(MPI_OK)

# --- reader
add_library(reader OBJECT reader.f90 reader_hdf5.f90 reader_ncdf.f90 reader_raw.f90)
target_link_libraries(reader PRIVATE const pathlib h5fortran::h5fortran nc4fortran::nc4fortran)

if(netcdf)
  add_executable(test_simsize_nc4 test_simsize_nc4.f90)
  target_link_libraries(test_simsize_nc4 nc4fortran::nc4fortran)
  # add_test(NAME unit:simsize_nc4 COMMAND test_simsize_nc4)
endif(netcdf)