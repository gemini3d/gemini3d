add_library(fsutils mkdir.F90 expanduser.f90 formats.f90)
target_compile_options(fsutils PRIVATE ${FFLAGS})
set_target_properties(fsutils PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

add_library(logging logging.f90)
target_compile_options(logging PRIVATE ${FFLAGS})
set_target_properties(logging PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
