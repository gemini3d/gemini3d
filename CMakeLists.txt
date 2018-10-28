cmake_minimum_required (VERSION 3.11)
project(gemini3d Fortran C)
enable_testing()

if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.12.0")  # quotes are needed!
  cmake_policy(SET CMP0074 NEW)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules/)

set(LIB_DIR ${PROJECT_SOURCE_DIR}/../fortran-libs) # for self-compiled libraries

include(${PROJECT_SOURCE_DIR}/cmake/compilers.cmake)

include(${PROJECT_SOURCE_DIR}/cmake/libraries.cmake)

# --- Build Gemini
include_directories(${CMAKE_CURRENT_BINARY_DIR}/numerical)

add_library(msis vendor/msis00/msis00_gfortran.f)
target_compile_options(msis PRIVATE -std=legacy -w -fno-implicit-none)

add_subdirectory(numerical)

add_library(io io/io.f90 io/expanduser.f90)
target_link_libraries(io PRIVATE const mpimod grid calculus)

add_library(temporal temporal/temporal.f90)
target_link_libraries(temporal PRIVATE const mpimod grid ${MPI_Fortran_LIBRARIES})
target_include_directories(temporal PRIVATE ${MPI_Fortran_INCLUDE_DIRS})

add_library(neutral neutral/neutral.f90)
target_link_libraries(neutral PUBLIC interp io msis temporal ${MPI_Fortran_LIBRARIES})
target_include_directories(neutral PUBLIC ${MPI_Fortran_INCLUDE_DIRS}) #PUBLIC for ionization with ifort

add_library(ionization ionization/ionization.f90)
target_link_libraries(ionization PRIVATE neutral)

add_library(collision collisions/collisions.f90)
target_link_libraries(collision PRIVATE const)

add_library(advec numerical/advection/advec_mpi.f90)
target_link_libraries(advec PRIVATE const mpimod grid)

add_library(sources sources/sources.f90)
target_link_libraries(sources PRIVATE collision grid calculus)

# --- MUMPS stuff
add_library(potential_mumps numerical/potential/potential_mumps.f90)
target_link_libraries(potential_mumps PRIVATE grid calculus ${MUMPS_LIBRARIES} ${MPI_Fortran_LIBRARIES})
target_include_directories(potential_mumps PRIVATE ${MUMPS_INCLUDE_DIRS} ${MPI_Fortran_INCLUDE_DIRS})

add_library(potential_comm_mumps numerical/potential/potential_comm_mumps.f90)
target_link_libraries(potential_comm_mumps PRIVATE potential_mumps potentialBCs collision calculus grid const)

add_library(potentialBCs numerical/potential/boundary_conditions/potentialBCs_mumps.f90)
target_link_libraries(potentialBCs PRIVATE grid calculus io const ${MUMPS_LIBRARIES} ${SCALAPACK_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_include_directories(potentialBCs PRIVATE ${MUMPS_INCLUDE_DIRS})

# --- other libs
add_library(precipBCs ionization/boundary_conditions/precipBCs_mod.f90)
target_link_libraries(precipBCs PRIVATE const grid io interp temporal)

add_library(multifluid multifluid/multifluid.f90)
target_link_libraries(multifluid PRIVATE precipBCs mpimod grid ionization sources advec diffusion)

# --- Main GEMINI executable
add_executable(gemini gemini.f90)
target_link_libraries(gemini PRIVATE mpimod io neutral precipBCs temporal multifluid 
                      potentialBCs potential_comm_mumps const)


#--------magnetic field calculation executable---------------
add_executable(magcalc magcalc.f90)
target_link_libraries(magcalc PRIVATE mpimod const grid io temporal)

# --- self-tests
include(${PROJECT_SOURCE_DIR}/cmake/test2D.cmake)
