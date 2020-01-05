# based on: https://github.com/Kitware/VTK/blob/master/CMake/FindNetCDF.cmake
# in general, NetCDF requires C compiler even if only using Fortran

function(netcdf_c)

pkg_check_modules(NC netcdf)  # C / CXX

find_path(NetCDF_INCLUDE_DIR
  NAMES netcdf.h
  HINTS ${NC_INCLUDE_DIRS}
  DOC "NetCDF include directories")

if(NOT NetCDF_INCLUDE_DIR)
  return()
endif()

find_library(NetCDF_C_LIBRARY
  NAMES netcdf
  HINTS ${NC_LIBRARY_DIRS} ${NC_LIBDIR}
  DOC "NetCDF C library")

if(NOT NetCDF_C_LIBRARY)
  return()
endif()

set(NetCDF_C_FOUND true PARENT_SCOPE)
set(NetCDF_INCLUDE_DIR ${NetCDF_INCLUDE_DIR} PARENT_SCOPE)
set(NetCDF_LIBRARY ${NetCDF_C_LIBRARY} PARENT_SCOPE)

endfunction(netcdf_c)


function(netcdf_fortran)

if(NOT CMAKE_Fortran_COMPILER)
  return()
endif()

pkg_check_modules(NC netcdf-fortran)
if(NOT NC_FOUND) # homebrew
  pkg_check_modules(NC netcdf)
endif()

find_path(NetCDF_Fortran_INCLUDE_DIR
  names netcdf.mod
  HINTS ${NC_INCLUDE_DIRS}
  DOC "NetCDF Fortran Include")

if(NOT NetCDF_Fortran_INCLUDE_DIR)
  return()
endif()

find_library(NetCDF_Fortran_LIBRARY
  NAMES netcdff
  HINTS ${NC_LIBRARY_DIRS} ${NC_LIBDIR}
  DOC "NetCDF Fortran library")

if(NOT NetCDF_Fortran_LIBRARY)
  return()
endif()

set(NetCDF_Fortran_FOUND true PARENT_SCOPE)
set(NetCDF_INCLUDE_DIR ${NetCDF_INCLUDE_DIR} ${NetCDF_Fortran_INCLUDE_DIR} PARENT_SCOPE)
set(NetCDF_LIBRARY ${NetCDF_LIBRARY} ${NetCDF_Fortran_LIBRARY} PARENT_SCOPE)

endfunction(netcdf_fortran)

#============================================================
cmake_policy(VERSION 3.3)

find_package(PkgConfig)

netcdf_c()

if(NetCDF_LIBRARY AND Fortran IN_LIST NetCDF_FIND_COMPONENTS)
  netcdf_fortran()
endif()
mark_as_advanced(NetCDF_LIBRARY NetCDF_INCLUDE_DIR)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NetCDF
  REQUIRED_VARS NetCDF_LIBRARY NetCDF_INCLUDE_DIR
  HANDLE_COMPONENTS)

if(NetCDF_FOUND)
  set(NetCDF_INCLUDE_DIRS ${NetCDF_INCLUDE_DIR})
  set(NetCDF_LIBRARIES ${NetCDF_LIBRARY})
endif()
