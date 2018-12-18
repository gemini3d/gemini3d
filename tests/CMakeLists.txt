cmake_minimum_required (VERSION 3.12)  
# 3.9: MPI::MPI_Fortran, 3.11: Lapack,  3.12 find_ROOT
project(test_gem Fortran)
enable_testing()

get_directory_property(hasParent PARENT_DIRECTORY)
if(NOT hasParent)
  if(NOT realbits)
    set(realbits 64)
  endif()

  list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake/Modules/)
  include(../cmake/compilers.cmake)
  include(../cmake/libraries.cmake)
  include(${PROJECT_SOURCE_DIR}/../cmake/mumps.cmake)
endif()

find_package(MPI REQUIRED COMPONENTS Fortran)

add_executable(testmpi test_mpi.f90)
target_link_libraries(testmpi PRIVATE MPI::MPI_Fortran)
target_compile_options(testmpi PRIVATE ${FFLAGS})
add_test(NAME MPIexist COMMAND testmpi)
  
set_tests_properties(MPIexist PROPERTIES TIMEOUT 5)
#---------------


add_executable(testmumps test_mumps.f90)
target_include_directories(testmumps PRIVATE ${MUMPS_INCLUDE_DIRS})
target_compile_options(testmumps PRIVATE ${FFLAGS})
target_link_libraries(testmumps PRIVATE
                     ${MUMPS_LIBRARIES}
                     MPI::MPI_Fortran)
add_test(NAME MUMPS
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${CMAKE_CURRENT_BINARY_DIR}/testmumps
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  
set_tests_properties(MUMPS PROPERTIES TIMEOUT 5)
