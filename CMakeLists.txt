cmake_minimum_required (VERSION 3.7)
project(gemini3d Fortran C)
enable_testing()

# Test data directories
set(TEST_REFERENCE_DIR $ENV{HOME}/zettergmdata/testdata/2Dtest_files/2Dtest_output)
set(TEST_OUTDIR ${CMAKE_CURRENT_BINARY_DIR}/test2d)


if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.12.0")  # quotes are needed!
  cmake_policy(SET CMP0074 NEW)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules/)

set(LIB_DIR ${PROJECT_SOURCE_DIR}/../fortran-libs) # for self-compiled libraries (future: git submodule of gemini3d)

## ----- release or debug?

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  add_compile_options(-g -O0)
else()
  add_compile_options(-g -O3)
endif()

## ------- Compiler options
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
  # -r8  after literals are fixed to "e" or "wp"
  if(CMAKE_BUILD_TYPE STREQUAL Debug)
     add_compile_options(-debug extended -check all -heap-arrays)
  endif()
  add_compile_options(-warn nounused -traceback -fp-stack-check)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
  # -fdefault-real-8  after literals are fixed to "e" or "wp"
  add_compile_options(-mtune=native -fimplicit-none -Wall -Wpedantic -Wextra -fexceptions -fbacktrace -Wno-unused-dummy-argument -Wno-unused-variable -Wno-unused-function -I/opt/local/lib/)
# -fstack-protector
# -ffpe-trap=invalid,zero,overflow)#,underflow)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL PGI)
  add_compile_options(-Mallocatable=03)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Flang) 
  add_compile_options(-Mallocatable=03)
  link_libraries(-static-flang-libs)
endif()
# -------  Libraries-----------
# MKL
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
find_package(MKL COMPONENTS MPI)
set(MKLROOT $ENV{MKLROOT})
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
  include_directories(${MKL_INCLUDE_DIRS})
  set(BLA_VENDOR Intel)
else()
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
endif()

# MPI
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
  if(NOT MPI_Fortran_FOUND)
    find_package(MPI REQUIRED)
    add_compile_options(${MPI_Fortran_COMPILE_OPTIONS})
  endif()
else()
  find_package(MPI REQUIRED)
  add_compile_options(${MPI_Fortran_COMPILE_OPTIONS})
endif()
message(STATUS "MPI:  include: " ${MPI_Fortran_INCLUDE_DIRS} " lib: " ${MPI_Fortran_LIBRARIES} " MPI flags: " ${MPI_Fortran_COMPILE_OPTIONS} " compiler: " ${MPI_Fortran_COMPILER})


# SCALAPACK
list(APPEND SCALAPACK_ROOT ${LIB_DIR}/scalapack)
find_package(SCALAPACK REQUIRED)

# MUMPS -- use cmake -DMUMPS_ROOT= for hint.
# Intel MKL-compiled MUMPS requires at the linker for the main executable:
# mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 mkl_intel_lp64 mkl_intel_thread mkl_core
list(APPEND METIS_ROOT ${LIB_DIR}/metis)
list(APPEND Scotch_ROOT ${LIB_DIR}/scotch)
list(APPEND MUMPS_ROOT ${LIB_DIR}/MUMPS)
find_package(MUMPS REQUIRED COMPONENTS d)  # apt install libmumps-dev
#------ executable --------
include_directories(${CMAKE_CURRENT_BINARY_DIR}/numerical)

add_library(msis neutral/msis00/msis00_gfortran.f)
target_compile_options(msis PRIVATE -w -fno-implicit-none)

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

# ---
add_library(precipBCs ionization/boundary_conditions/precipBCs_mod.f90)
target_link_libraries(precipBCs PRIVATE const grid io interp temporal)

add_library(multifluid multifluid/multifluid.f90)
target_link_libraries(multifluid PRIVATE precipBCs mpimod grid ionization sources advec diffusion)

#------- main executable ------------
add_executable(gemini gemini.f90)

target_link_libraries(gemini PRIVATE mpimod io neutral precipBCs temporal multifluid 
                      potentialBCs potential_comm_mumps const ${MPI_Fortran_LIBRARIES})
target_include_directories(gemini PRIVATE ${MUMPS_INCLUDE_DIRS} ${MPI_Fortran_INCLUDE_DIRS})


#--------magnetic field calculation executable---------------
add_executable(magcalc magcalc.f90)
target_link_libraries(magcalc mpimod const grid io temporal)


# -- test main exe ---

add_test(NAME Gemini2D 
         COMMAND mpirun -np 4 ./gemini ${CMAKE_SOURCE_DIR}/initialize/2Dtest/config.ini ${TEST_OUTDIR} 
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
set_tests_properties(Gemini2D PROPERTIES 
                     TIMEOUT 900  # test should complete in under 15 minutes on ~ 2014 dual-core laptop
                     FIXTURES_SETUP GemMain)
                     

# ----- evaluate output accuracy vs. reference from Matt's HPC
find_package(Octave)
if (OCTAVE_MAJOR_VERSION GREATER_EQUAL 4)

add_test(NAME ReferenceCompare 
         COMMAND octave-cli -q --eval "compare_all('${TEST_OUTDIR}','${TEST_REFERENCE_DIR}')"
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests)
set_tests_properties(ReferenceCompare PROPERTIES
                     TIMEOUT 5)
                     #FIXTURES_REQUIRED GemMain)

endif()

# check that *.m syntax is Matlab compatible
find_package(Matlab COMPONENTS MAIN_PROGRAM)
if (Matlab_FOUND)
add_test(NAME MatlabCompare 
         COMMAND matlab -nojvm -r "exit(compare_all('${TEST_OUTDIR}','${TEST_REFERENCE_DIR}'))"
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests)
set_tests_properties(MatlabCompare PROPERTIES
                     TIMEOUT 30)  # Matlab with a lot of toolboxes takes ~ 15 seconds just to start
                     #FIXTURES_REQUIRED GemMain)

endif()
