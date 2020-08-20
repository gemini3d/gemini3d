set(CDEFS "Add_")
# "Add_" works for all modern compilers we tried.

add_library(mpi_stub elapse.c mpi.f mpic.c)
target_compile_definitions(mpi_stub PRIVATE ${CDEFS})
target_include_directories(mpi_stub INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(mpi_stub PROPERTIES
  Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}
  ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  target_compile_options(mpi_stub PRIVATE -std=legacy)
endif()
