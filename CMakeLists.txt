cmake_minimum_required(VERSION 3.12)  # 3.12 for find_ROOT

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Debug or Release" FORCE)
endif()

project(gemini3d
  LANGUAGES C Fortran  # MUST include C language for Intel / MKL to work
  HOMEPAGE_URL https://github.com/gemini3d/gemini)

enable_testing()

if(NOT realbits)
  set(realbits 64)
endif()

option(TRACE "dump variables to disk at certain locations in program" OFF)
option(useglow "use NCAR GLOW instead of Fang")
option(usematab "enable Matlab tests, which are slow and duplicate fast Octave tests" OFF)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules/)

include(${PROJECT_SOURCE_DIR}/cmake/compilers.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/libraries.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/mumps.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/utils.cmake)

# --- pretests, to help catch missing libs
add_subdirectory(tests)

# --- vendor libraries
add_subdirectory(vendor/msis00)

# --- GEMINI
add_subdirectory(ionization)

add_subdirectory(numerical)

add_subdirectory(io)
add_library(io io/io.F90)
target_link_libraries(io PRIVATE const mpimod grid calculus fsutils logging)
target_compile_definitions(io PRIVATE TRACE=${TRACE})
target_compile_options(io PRIVATE ${FFLAGS})

add_library(temporal temporal/temporal.f90)
target_link_libraries(temporal PRIVATE const mpimod grid MPI::MPI_Fortran)
target_compile_options(temporal PRIVATE ${FFLAGS})

add_subdirectory(temporal)

add_library(neutral neutral/neutral.f90)
target_link_libraries(neutral PUBLIC interp io msis timeutils MPI::MPI_Fortran)
target_compile_options(neutral PRIVATE ${FFLAGS})

add_library(collision collisions/collisions.f90)
target_link_libraries(collision PRIVATE const)
target_compile_options(collision PRIVATE ${FFLAGS})

add_library(sources sources/sources.f90)
target_link_libraries(sources PRIVATE collision grid calculus)
target_compile_options(sources PRIVATE ${FFLAGS})

# --- MUMPS stuff
add_library(potential_mumps numerical/potential/potential_mumps.F90)
target_link_libraries(potential_mumps PRIVATE const grid calculus interp
  ${MUMPS_LIBRARIES} MPI::MPI_Fortran)
target_compile_definitions(potential_mumps PRIVATE REALBITS=${realbits})
target_include_directories(potential_mumps PRIVATE ${MUMPS_INCLUDE_DIRS})
target_compile_options(potential_mumps PRIVATE ${FFLAGS})

add_library(potential_comm_mumps numerical/potential/potential_comm_mumps.f90)
target_link_libraries(potential_comm_mumps PRIVATE potential_mumps potentialBCs collision calculus grid const)
target_compile_options(potential_comm_mumps PRIVATE ${FFLAGS})

add_library(potentialBCs numerical/potential/boundary_conditions/potentialBCs_mumps.f90)
target_link_libraries(potentialBCs PRIVATE temporal timeutils grid calculus interp io const
    ${MUMPS_LIBRARIES} MPI::MPI_Fortran)
target_include_directories(potentialBCs PRIVATE ${MUMPS_INCLUDE_DIRS})
target_compile_options(potentialBCs PRIVATE ${FFLAGS})

# --- other libs
add_library(multifluid multifluid/multifluid.f90)
target_link_libraries(multifluid PRIVATE precipBCs mpimod grid ionization sources collision advec diffusion)
target_compile_options(multifluid PRIVATE ${FFLAGS})

# --- Main GEMINI executable
add_executable(gemini_fang.bin gemini.f90)
target_link_libraries(gemini_fang.bin PRIVATE mpimod io neutral precipBCs
  temporal multifluid
  potentialBCs potential_comm_mumps const
  ${LAPACK_LIBRARIES})  # necessary for self-compiled
#  MPI::MPI_Fortran)
target_compile_options(gemini_fang.bin PRIVATE ${FFLAGS})

if(useglow)
  add_library(multifluid_glow multifluid/multifluid.f90)
  target_link_libraries(multifluid_glow PRIVATE precipBCs mpimod grid ionization sources collision advec diffusion)
  target_compile_options(multifluid_glow PRIVATE ${FFLAGS})

  add_executable(gemini_glow.bin gemini.f90)
  target_link_libraries(gemini_glow.bin PRIVATE mpimod io neutral precipBCs
    temporal multifluid_glow
    potentialBCs potential_comm_mumps const
    ${LAPACK_LIBRARIES})  # necessary for self-compiled
  #  MPI::MPI_Fortran)
  target_compile_options(gemini_glow.bin PRIVATE ${FFLAGS})
endif(useglow)

#--------magnetic field calculation executable---------------
add_executable(magcalc.bin magcalc.f90)
target_link_libraries(magcalc.bin PRIVATE mpimod const grid io timeutils)
target_compile_options(magcalc.bin PRIVATE ${FFLAGS})

# --- self-tests
include(${PROJECT_SOURCE_DIR}/cmake/test2d.cmake)
if(NOT DEFINED ENV{CI} OR NOT $ENV{CI})
  # need the dollar sign on ENV for IF: save needless downloading
  include(${PROJECT_SOURCE_DIR}/cmake/test3d.cmake)
endif()
