add_library(fsutils mkdir.F90 expanduser.f90)

add_library(logging logging.f90)

#=== io library
if(hdf5)
  set(io_src plasma_output_hdf5.f90 mag_hdf5.f90 aurora_hdf5.f90)
elseif(netcdf)
  set(io_src plasma_output_ncdf.f90 mag_ncdf.f90 aurora_ncdf.f90)
else()
  set(io_src plasma_output_raw.f90 mag_raw.f90 aurora_raw.f90)
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  list(APPEND io_src path_exists_intel.f90)
else()
  list(APPEND io_src path_exists.f90)
endif()

add_library(io io.f90 aurora.f90 output.f90 input.f90
  mag.f90
  plasma.f90 plasma_input_raw.f90
  ${io_src})
target_link_libraries(io PRIVATE calculus const fsutils grid mpimod timeutils)

if(hdf5)
  target_link_libraries(io PUBLIC h5fortran::h5fortran)
elseif(netcdf)
  target_include_directories(io PUBLIC ${NetCDF_INCLUDE_DIRS})
  target_link_libraries(io PUBLIC nc4fortran ${NetCDF_LIBRARIES})
endif()
