add_library(logging OBJECT logging.f90)

# --- pathlib
add_library(pathlib OBJECT pathlib.f90)
if(WIN32)
  target_sources(pathlib PRIVATE pathlib_windows.f90)
else()
  target_sources(pathlib PRIVATE pathlib_unix.f90)
endif()
if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  target_sources(pathlib PRIVATE path_exists_intel.f90)
else()
  target_sources(pathlib PRIVATE path_exists.f90)
endif()

# --- config
add_library(config OBJECT config.f90)
target_link_libraries(config PRIVATE pathlib const)

add_executable(test_config test_config.f90)
target_link_libraries(test_config PRIVATE config const pathlib)
add_test(NAME unit::config COMMAND $<TARGET_FILE:test_config> ${PROJECT_SOURCE_DIR}/tests/test2deq_config.nml)
set_tests_properties(unit::config PROPERTIES
REQUIRED_FILES ${PROJECT_SOURCE_DIR}/tests/test2deq_config.nml
TIMEOUT 10)


# --- io
add_library(io OBJECT io.f90 output.f90
aurora.f90 aurora_raw.f90 aurora_hdf5.f90 aurora_ncdf.f90
mag.f90
plasma.f90)

target_link_libraries(io PRIVATE calculus const reader grid mpimod timeutils)

if(hdf5)
  target_sources(io PRIVATE plasma_input_hdf5.f90 plasma_output_hdf5.f90 mag_hdf5.f90)
  target_link_libraries(io PUBLIC h5fortran::h5fortran)
endif()

if(netcdf)
  target_sources(io PRIVATE plasma_input_ncdf.f90 plasma_output_ncdf.f90 mag_ncdf.f90)
  target_include_directories(io PUBLIC ${NetCDF_INCLUDE_DIRS})
  target_link_libraries(io PUBLIC nc4fortran ${NetCDF_LIBRARIES})
endif()

if(NOT hdf5 AND NOT netcdf)
  target_sources(io PRIVATE plasma_input_raw.f90 plasma_output_raw.f90 mag_raw.f90)
endif()


# --- reader
add_library(reader OBJECT reader.f90)
target_link_libraries(reader PRIVATE const)

if(hdf5)
  target_sources(reader PRIVATE reader_hdf5.f90)
  target_link_libraries(reader PRIVATE h5fortran::h5fortran)
endif()

if(netcdf)
  target_sources(reader PRIVATE reader_nc4.f90)
  target_link_libraries(reader PRIVATE nc4fortran::nc4fortran)
endif()

if(NOT hdf5 AND NOT netcdf)
  target_sources(reader PRIVATE reader_raw.f90)
endif()
