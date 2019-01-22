cmake_minimum_required(VERSION 3.12)
# 3.9: MPI::MPI_Fortran  3.11 for FindLapack bug   3.12 for find_ROOT
project(gemini3d C Fortran)  # MUST include C language for Intel / MKL to work
enable_testing()

if(NOT realbits)
  set(realbits 64)
endif()

option(TRACE "dump variables to disk at certain locations in program" off)
option(USEGLOW "use NCAR GLOW instead of Fang")

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
add_library(ionization ionization/ionization.f90)
target_link_libraries(ionization PRIVATE neutral const)
target_include_directories(ionization PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/numerical)
if(USEGLOW)
  include(${PROJECT_SOURCE_DIR}/cmake/glow.cmake)

  target_sources(ionization PRIVATE ionization/glow_run.f90)
  target_link_libraries(ionization PRIVATE cglow)
  target_include_directories(ionization PRIVATE ${ncarglow_BINARY_DIR})
else()
  target_sources(ionization PRIVATE ionization/glow_dummy.f90)
endif(USEGLOW)

add_subdirectory(numerical)

add_library(io io/io.F90 io/expanduser.f90)
target_link_libraries(io PRIVATE const mpimod grid calculus)
target_compile_definitions(io PRIVATE TRACE=${TRACE})
target_include_directories(io PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/numerical)
target_compile_options(io PRIVATE ${FFLAGS})

add_library(temporal temporal/temporal.f90)
target_link_libraries(temporal PRIVATE const mpimod grid MPI::MPI_Fortran)
target_include_directories(temporal PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/numerical)
target_compile_options(temporal PRIVATE ${FFLAGS})

add_library(timeutils temporal/timeutils.f90)
target_link_libraries(timeutils PRIVATE const)
target_include_directories(timeutils PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/numerical)
target_compile_options(timeutils PRIVATE ${FFLAGS})

add_library(neutral neutral/neutral.f90)
target_link_libraries(neutral PUBLIC interp io msis timeutils MPI::MPI_Fortran)
target_include_directories(neutral PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/numerical)
target_compile_options(neutral PRIVATE ${FFLAGS})



add_library(collision collisions/collisions.f90)
target_link_libraries(collision PRIVATE const)
target_include_directories(collision PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/numerical)
target_compile_options(collision PRIVATE ${FFLAGS})

add_library(sources sources/sources.f90)
target_link_libraries(sources PRIVATE collision grid calculus)
target_include_directories(sources PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/numerical)
target_compile_options(sources PRIVATE ${FFLAGS})

# --- MUMPS stuff
add_library(potential_mumps numerical/potential/potential_mumps.F90)
target_link_libraries(potential_mumps PRIVATE const grid calculus
       ${MUMPS_LIBRARIES}
       MPI::MPI_Fortran)
target_compile_definitions(potential_mumps PRIVATE REALBITS=${realbits})
target_include_directories(potential_mumps PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/numerical ${MUMPS_INCLUDE_DIRS})
target_compile_options(potential_mumps PRIVATE ${FFLAGS})

add_library(potential_comm_mumps numerical/potential/potential_comm_mumps.f90)
target_link_libraries(potential_comm_mumps PRIVATE potential_mumps potentialBCs collision calculus grid const)
target_include_directories(potential_comm_mumps PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/numerical)
target_compile_options(potential_comm_mumps PRIVATE ${FFLAGS})

add_library(potentialBCs numerical/potential/boundary_conditions/potentialBCs_mumps.f90)
target_link_libraries(potentialBCs PRIVATE temporal grid calculus interp io const
    ${MUMPS_LIBRARIES}
#    ${BLAS_LIBRARIES}
    MPI::MPI_Fortran)
target_include_directories(potentialBCs PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/numerical
  ${MUMPS_INCLUDE_DIRS})
  target_compile_options(potentialBCs PRIVATE ${FFLAGS})

# --- other libs
add_library(precipBCs ionization/boundary_conditions/precipBCs_mod.f90)
target_link_libraries(precipBCs PRIVATE const grid io interp temporal)
target_include_directories(precipBCs PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/numerical)
target_compile_options(precipBCs PRIVATE ${FFLAGS})

add_library(multifluid multifluid/multifluid.f90)
target_link_libraries(multifluid PRIVATE precipBCs mpimod grid ionization sources collision advec diffusion)
target_include_directories(multifluid PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/numerical)
target_compile_options(multifluid PRIVATE ${FFLAGS})

# --- Main GEMINI executable
add_executable(gemini.bin gemini.f90)
target_link_libraries(gemini.bin PRIVATE mpimod io neutral precipBCs
  temporal multifluid
  potentialBCs potential_comm_mumps const)
#  MPI::MPI_Fortran)
target_include_directories(gemini.bin PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/numerical)
target_compile_options(gemini.bin PRIVATE ${FFLAGS})


#--------magnetic field calculation executable---------------
add_executable(magcalc.bin magcalc.f90)
target_link_libraries(magcalc.bin PRIVATE mpimod const grid io timeutils)
target_include_directories(magcalc.bin PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/numerical)
target_compile_options(magcalc.bin PRIVATE ${FFLAGS})

# --- self-tests
if(USEGLOW)
  include(${PROJECT_SOURCE_DIR}/cmake/test2d_glow.cmake)
  include(${PROJECT_SOURCE_DIR}/cmake/test3d_glow.cmake)
else()
  include(${PROJECT_SOURCE_DIR}/cmake/test2d.cmake)
  include(${PROJECT_SOURCE_DIR}/cmake/test3d.cmake)
endif()
