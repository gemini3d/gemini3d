find_package(NetCDF REQUIRED COMPONENTS Fortran)

# -- verify

set(CMAKE_REQUIRED_INCLUDES ${NetCDF_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${NetCDF_LIBRARIES})

include(CheckFortranSourceCompiles)
check_fortran_source_compiles("use netcdf, only : nf90_open; end" NetCDF_OK SRC_EXT f90)

if(NOT NetCDF_OK)
  message(FATAL_ERROR "NetCDF library not working with ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
endif(NOT NetCDF_OK)
