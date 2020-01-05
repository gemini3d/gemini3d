cmake_minimum_required(VERSION 3.14)  # 3.12 for find_ROOT, 3.13 for a lot of user functions, 3.14 for check_fortran_source_runs

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Debug or Release" FORCE)
endif()

project(gemini3d
  LANGUAGES C Fortran  # MUST include C language for Intel / MKL to work
  HOMEPAGE_URL https://github.com/gemini3d/gemini)

enable_testing()

if(WIN32 AND NOT CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  message(FATAL_ERROR "On native Windows, only the Intel compilers (ifort) works at this time, because OpenMPI is not available for Windows. Try Windows Subsystem for Linux.")
endif()

if(NOT realbits)
  set(realbits 64)
endif()

option(glow "use NCAR GLOW instead of Fang" ON)
option(matlab "enable Matlab tests, which are slow and duplicate fast Octave tests" OFF)
option(hdf5 "use HDF5 instead of raw uncompressed output")
option(metis "use METIS" OFF)
option(scotch "use Scotch" OFF)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules/)

include(${CMAKE_SOURCE_DIR}/cmake/compilers.cmake)
if(hdf5)
include(${CMAKE_SOURCE_DIR}/cmake/h5fortran.cmake)
endif()
include(${CMAKE_SOURCE_DIR}/cmake/utils.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/mumps.cmake)

# --- pretests, to help catch missing libs
add_subdirectory(tests)

# --- vendor libraries
add_subdirectory(vendor/msis00)

# --- GEMINI
add_subdirectory(ionization)

add_subdirectory(numerical)

add_subdirectory(io)

add_subdirectory(temporal)

add_library(neutral neutral/neutral.f90)
target_link_libraries(neutral PRIVATE const grid interp mpimod msis timeutils MPI::MPI_Fortran)
target_compile_options(neutral PRIVATE ${FFLAGS})

add_library(collision collisions/collisions.f90)
target_link_libraries(collision PRIVATE const)
target_compile_options(collision PRIVATE ${FFLAGS})

add_library(sources sources/sources.f90)
target_link_libraries(sources PRIVATE calculus collision const grid mpimod)
target_compile_options(sources PRIVATE ${FFLAGS})

# --- other libs
add_library(multifluid multifluid/multifluid.f90)
target_link_libraries(multifluid PRIVATE advec calculus collision const diffusion grid ionization mpimod precipBCs sources timeutils)
target_compile_options(multifluid PRIVATE ${FFLAGS})

# --- Main GEMINI executable
add_executable(gemini_fang.bin gemini.f90)
target_link_libraries(gemini_fang.bin PRIVATE const grid io mpimod multifluid neutral potential precipBCs temporal timeutils
  ${LAPACK_LIBRARIES})  # necessary for self-compiled
target_compile_options(gemini_fang.bin PRIVATE ${FFLAGS})

if(glow)
  add_library(multifluid_glow multifluid/multifluid.f90)
  target_link_libraries(multifluid_glow PRIVATE advec calculus collision const diffusion grid ionization_glow mpimod precipBCs sources timeutils)
  target_compile_options(multifluid_glow PRIVATE ${FFLAGS})

  add_executable(gemini_glow.bin gemini.f90)
  target_link_libraries(gemini_glow.bin PRIVATE const grid io mpimod multifluid_glow neutral potential precipBCs temporal timeutils
    ${LAPACK_LIBRARIES})  # necessary for self-compiled
  target_compile_options(gemini_glow.bin PRIVATE ${FFLAGS})
  target_link_directories(gemini_glow.bin PRIVATE ${GLOW_BINARY_DIR})  # "not found -lcglow" without this
endif(glow)

#--------magnetic field calculation executable---------------
add_executable(magcalc.bin magcalc.f90)
target_link_libraries(magcalc.bin PRIVATE mpimod const grid io timeutils)
target_compile_options(magcalc.bin PRIVATE ${FFLAGS})

# --- self-tests
include(${CMAKE_SOURCE_DIR}/cmake/test2d.cmake)
#if(NOT DEFINED ENV{CI} OR NOT $ENV{CI})
  # need the dollar sign on ENV for IF: save needless downloading
include(${CMAKE_SOURCE_DIR}/cmake/test3d.cmake)
#endif()
