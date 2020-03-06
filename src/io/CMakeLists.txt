add_library(pathlib OBJECT pathlib.f90)
if(WIN32)
  target_sources(pathlib PRIVATE pathlib_windows.f90)
else()
  target_sources(pathlib PRIVATE pathlib_unix.f90)
endif()

add_library(logging OBJECT logging.f90)

#=== io library
add_library(io OBJECT io.f90 output.f90 input.f90
aurora.f90
mag.f90
plasma.f90)

target_link_libraries(io PRIVATE calculus const pathlib reader grid mpimod timeutils)

if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  target_sources(io PRIVATE path_exists_intel.f90)
else()
  target_sources(io PRIVATE path_exists.f90)
endif()


if(hdf5)
  target_sources(io PRIVATE plasma_input_hdf5.f90 plasma_output_hdf5.f90 mag_hdf5.f90 aurora_hdf5.f90)
endif()

if(netcdf)
  target_sources(io PRIVATE plasma_input_ncdf.f90 plasma_output_ncdf.f90 mag_ncdf.f90 aurora_ncdf.f90)
endif()

if(NOT hdf5 AND NOT netcdf)
  target_sources(io PRIVATE plasma_input_raw.f90 plasma_output_raw.f90 mag_raw.f90 aurora_raw.f90)
endif()


if(hdf5)
  target_link_libraries(io PUBLIC h5fortran::h5fortran)
elseif(netcdf)
  target_include_directories(io PUBLIC ${NetCDF_INCLUDE_DIRS})
  target_link_libraries(io PUBLIC nc4fortran ${NetCDF_LIBRARIES})
endif()

#=== reader
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
