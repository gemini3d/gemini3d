add_library(pathlib pathlib.F90)
target_include_directories(pathlib INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(pathlib PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

add_library(logging logging.f90)

#=== io library
if(hdf5)
  set(io_src plasma_input_hdf5.f90 plasma_output_hdf5.f90 mag_hdf5.f90 aurora_hdf5.f90)
elseif(netcdf)
  set(io_src plasma_input_ncdf.f90 plasma_output_ncdf.f90 mag_ncdf.f90 aurora_ncdf.f90)
else()
  set(io_src plasma_input_raw.f90 plasma_output_raw.f90 mag_raw.f90 aurora_raw.f90)
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  list(APPEND io_src path_exists_intel.f90)
else()
  list(APPEND io_src path_exists.f90)
endif()

add_library(io io.f90 aurora.f90 output.f90 input.f90 mag.f90 plasma.f90 ${io_src})
target_link_libraries(io PRIVATE calculus const pathlib reader grid mpimod timeutils)
target_include_directories(io INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(io PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

if(hdf5)
  target_link_libraries(io PUBLIC h5fortran::h5fortran)
elseif(netcdf)
  target_include_directories(io PUBLIC ${NetCDF_INCLUDE_DIRS})
  target_link_libraries(io PUBLIC nc4fortran ${NetCDF_LIBRARIES})
endif()

#=== reader
add_library(reader reader.f90)
if(hdf5)
  target_sources(reader PRIVATE reader_hdf5.f90)
  target_link_libraries(reader PRIVATE h5fortran::h5fortran)
elseif(netcdf)
  target_sources(reader PRIVATE reader_nc4.f90)
  target_link_libraries(reader PRIVATE nc4fortran::nc4fortran)
else()
  target_sources(reader PRIVATE reader_raw.f90)
endif()
target_link_libraries(reader PRIVATE const)
target_include_directories(reader INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(reader PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)