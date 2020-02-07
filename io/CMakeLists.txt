add_library(pathlib pathlib.F90)
target_include_directories(pathlib INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(pathlib PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

# FIXME: test in Fortran not C
if(WIN32)
  set(_code "
#include <stdlib.h>
#include <string.h>
int main(void) {
char resolved_path[_MAX_PATH];
_fullpath(resolved_path, \"..\", _MAX_PATH);
if(strcmp(\"..\", resolved_path) == 0 || strcmp(\".\", resolved_path) == 0) return 1;
return 0; }")
  check_c_source_runs("${_code}" windowsio_OK)
  set(posixio_OK false)
elseif(UNIX AND NOT APPLE)
  set(_code "
#include <limits.h>
#include <stdlib.h>
#include <string.h>
int main(void) {
char resolved_path[PATH_MAX];
realpath(\"..\", resolved_path);
if(strcmp(\"..\", resolved_path) == 0 || strcmp(\".\", resolved_path) == 0) return 1;
return 0; }")
  check_c_source_runs("${_code}" posixio_OK)
  set(windowsio_OK false)
else()
  set(windowsio_OK false)
  set(posixio_OK false)
endif()

if(posixio_OK)
  target_sources(pathlib PRIVATE io_posix.f90)
elseif(windowsio_OK)
  target_sources(pathlib PRIVATE io_windows.f90)
else()
  target_sources(pathlib PRIVATE io_dummy.f90)
endif()

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
target_link_libraries(io PRIVATE calculus const pathlib grid mpimod timeutils)
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