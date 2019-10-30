add_library(fsutils mkdir.F90 expanduser.f90)
target_compile_options(fsutils PRIVATE ${FFLAGS})
set_target_properties(fsutils PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

add_library(logging logging.f90)
target_compile_options(logging PRIVATE ${FFLAGS})
set_target_properties(logging PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

if(usehdf5)
  set(io_src plasma_output_hdf5.f90 mag_hdf5.f90)
  add_dependencies(io h5fortran)
  target_link_directories(io PRIVATE ${h5fortran_BINARY_DIR})
  target_include_directories(io PRIVATE ${h5fortran_BINARY_DIR})
else()
  set(io_src plasma_output_raw.f90 mag_raw.f90)
endif()
add_library(io io.f90 aurora.f90 output.f90 input.f90
  mag.f90
  plasma.f90 plasma_input_raw.f90
  ${io_src})
set_target_properties(io PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
target_link_libraries(io PRIVATE calculus const fsutils grid mpimod timeutils)
target_compile_definitions(io PRIVATE TRACE=${TRACE})
target_compile_options(io PRIVATE ${FFLAGS})