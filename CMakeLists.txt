cmake_minimum_required(VERSION 3.14)  # 3.12 for find_ROOT, 3.13 for a lot of user functions, 3.14 for check_fortran_source_runs & fetchcontent

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Debug or Release" FORCE)
endif()

project(gemini3d
  LANGUAGES C Fortran  # MUST include C language for Intel / MKL to work
  DESCRIPTION "3-D ionospheric model"
  HOMEPAGE_URL https://github.com/gemini3d/gemini)

enable_testing()

if(NOT realbits)
  set(realbits 64)
endif()
option(glow "use NCAR GLOW instead of Fang" ON)
option(matlab "enable Matlab tests, which are slow and duplicate fast Octave tests" OFF)
option(hdf5 "use HDF5 file I/O")
option(netcdf "use NetCDF file I/O" OFF)
option(metis "use METIS" OFF)
option(scotch "use Scotch" OFF)

if(netcdf AND hdf5)
  message(FATAL_ERROR "HDF5 and NetCDF are mutually exclusive")
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules/)

include(${CMAKE_SOURCE_DIR}/cmake/compilers.cmake)
if(hdf5)
include(${CMAKE_SOURCE_DIR}/cmake/h5fortran.cmake)
elseif(netcdf)
include(${CMAKE_SOURCE_DIR}/cmake/nc4fortran.cmake)
endif()
include(${CMAKE_SOURCE_DIR}/cmake/utils.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/mumps.cmake)

# --- simple tests to help catch broken libs
add_subdirectory(tests)

# --- vendor libraries
add_subdirectory(vendor/msis00)

# --- GEMINI
add_subdirectory(ionization)
add_subdirectory(numerical)
add_subdirectory(io)
add_subdirectory(temporal)
add_subdirectory(neutral)

add_library(collision collisions/collisions.f90)
target_link_libraries(collision PRIVATE const)
target_include_directories(collision INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(collision PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

add_subdirectory(sources)

# --- other libs
add_library(multifluid multifluid/multifluid.f90)
target_link_libraries(multifluid PRIVATE advec calculus collision const diffusion grid mesh ionization mpimod precipBCs sources timeutils)

# --- Main GEMINI executable
add_executable(gemini.bin gemini.f90)
target_link_libraries(gemini.bin PRIVATE const grid mesh io mpimod multifluid neutral potential precipBCs temporal timeutils
  ${LAPACK_LIBRARIES})  # necessary for self-compiled

#--------magnetic field calculation executable---------------
add_executable(magcalc.bin magcalc.f90)
target_link_libraries(magcalc.bin PRIVATE mpimod const grid mesh io timeutils)

# --- self-tests
include(${CMAKE_SOURCE_DIR}/cmake/test2d.cmake)
#if(NOT DEFINED ENV{CI} OR NOT $ENV{CI})
  # need the dollar sign on ENV for IF: save needless downloading
include(${CMAKE_SOURCE_DIR}/cmake/test3d.cmake)
#endif()
