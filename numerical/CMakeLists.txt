cmake_minimum_required (VERSION 3.11)
project(gem-num Fortran)
enable_testing()

if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.12.0")  # quotes are needed!
  cmake_policy(SET CMP0074 NEW)
endif()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake/Modules/)
## Lapack95
list(APPEND LAPACK95_ROOT ${PROJECT_SOURCE_DIR}/../../fortran-libs/LAPACK95)
find_package(LAPACK95)
if(NOT LAPACK95_FOUND)
  find_package(LAPACK REQUIRED)
  add_library(f95_lapack ${PROJECT_SOURCE_DIR}/../vendor/lapack95/gbsv.f90)
  target_link_libraries(f95_lapack PRIVATE const)
  set(LAPACK95_LIBRARIES f95_lapack ${LAPACK_LIBRARIES})
endif()

#include_directories(${CMAKE_SOURCE_DIR}/objects)

add_library(const constants/phys_consts.f90)
#set_target_properties(const PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_SOURCE_DIR}/objects)

add_subdirectory(interpolation)

add_library(mpimod mpimod/mpimod.f90)
target_link_libraries(mpimod PUBLIC const ${MPI_Fortran_LIBRARIES})
target_include_directories(mpimod PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
#set_target_properties(mpimod PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_SOURCE_DIR}/objects)


add_library(grid grid/grid.f90)
target_link_libraries(grid PUBLIC const mpimod ${MPI_Fortran_LIBRARIES})
target_include_directories(grid PUBLIC ${MPI_Fortran_INCLUDE_DIRS})  # PUBLIC for calculus on Intel 
#set_target_properties(grid PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_SOURCE_DIR}/objects)


add_library(calculus calculus/calculus.f90)
target_link_libraries(calculus PRIVATE grid)
#set_target_properties(calculus PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_SOURCE_DIR}/objects)


# ----- diffusion
add_library(diffusion)  # CMake >= 3.11
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
  target_sources(diffusion PRIVATE diffusion/diffusion_MKL.f90)
else()
  target_sources(diffusion PRIVATE diffusion/diffusion.f90)
endif()

target_link_libraries(diffusion PRIVATE grid ${LAPACK95_LIBRARIES})
target_include_directories(diffusion PRIVATE ${LAPACK95_INCLUDE_DIRS})
#set_target_properties(diffusion PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_SOURCE_DIR}/objects)


add_executable(test_diffusion1D diffusion/test_diffusion1D.f90)
target_link_libraries(test_diffusion1D diffusion const)

add_test(NAME diffusion1D COMMAND test_diffusion1D)
set_tests_properties(diffusion1D PROPERTIES
                     TIMEOUT 5
                     FIXTURES_SETUP GemDiff)

find_package(Octave)
if (OCTAVE_MAJOR_VERSION GREATER_EQUAL 4)

add_test(NAME OctaveDiffusion1D 
         COMMAND octave-cli -q --eval "test_diffusion1D('${CMAKE_CURRENT_BINARY_DIR}/test_diffusion1d.dat')"
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/diffusion)

set_tests_properties(OctaveDiffusion1D  PROPERTIES
                     TIMEOUT 5
                     FIXTURES_REQUIRED GemDiff)
          
endif()

