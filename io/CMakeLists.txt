add_library(fsutils mkdir.F90 expanduser.f90)
target_compile_options(fsutils PRIVATE ${FFLAGS})
set_target_properties(fsutils PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

add_library(logging logging.f90)
target_compile_options(logging PRIVATE ${FFLAGS})
set_target_properties(logging PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

#=== io library
if(hdf5)
  set(io_src plasma_output_hdf5.f90 mag_hdf5.f90 aurora_hdf5.f90)
else()
  set(io_src plasma_output_raw.f90 mag_raw.f90 aurora_raw.f90)
endif()

add_library(io io.f90 aurora.f90 output.f90 input.F90
  mag.f90
  plasma.f90 plasma_input_raw.f90
  ${io_src})
set_target_properties(io PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
target_link_libraries(io PRIVATE calculus const fsutils grid mpimod timeutils)
target_compile_definitions(io PRIVATE TRACE=${TRACE})
target_compile_options(io PRIVATE ${FFLAGS})

if(hdf5)
  # these must be PUBLIC so that the top-level .bin link properly
  target_link_directories(io PUBLIC ${h5fortran_BINARY_DIR})
  target_include_directories(io PUBLIC ${h5fortran_BINARY_DIR} ${HDF5_INCLUDE_DIRS} ${HDF5_Fortran_INCLUDE_DIRS})
  target_compile_options(io PUBLIC ${HDF5_Fortran_DEFINITIONS})
  add_dependencies(io h5fortran_proj)
  target_link_libraries(io PUBLIC h5fortran ${HDF5_Fortran_HL_LIBRARIES} ${HDF5_Fortran_LIBRARIES})
endif()
