cmake_minimum_required (VERSION 3.0)
project(test_gem Fortran C)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake/Modules/)
enable_testing()

 
#------------ MPI

find_package(MPI REQUIRED COMPONENTS Fortran)

find_package(Threads REQUIRED)
add_compile_options(${MPI_Fortran_COMPILE_OPTIONS})
include_directories(${MPI_Fortran_INCLUDE_DIRS})

add_executable(testmpi test_mpi.f90)
target_link_libraries(testmpi PRIVATE ${MPI_Fortran_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

add_test(NAME MPIexist COMMAND testmpi)

#------------- 
