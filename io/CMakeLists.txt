cmake_minimum_required(VERSION 3.1)
project(filesystemutils Fortran)
enable_testing()

add_library(fsutils mkdir.F90 expanduser.f90 formats.f90)
target_compile_options(fsutils PRIVATE ${FFLAGS})
set_target_properties(fsutils PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(test_formats test_formats.f90)
target_compile_options(test_formats PRIVATE ${FFLAGS})
target_link_libraries(test_formats PRIVATE fsutils)
set_target_properties(test_formats PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME DateFormats COMMAND test_formats)