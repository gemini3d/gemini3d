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

if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  list(APPEND io_src path_exists_intel.f90)
else()
  list(APPEND io_src path_exists.f90)
endif()

add_library(io io.f90 aurora.f90 output.f90 input.f90
  mag.f90
  plasma.f90 plasma_input_raw.f90
  ${io_src})
set_target_properties(io PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
target_link_libraries(io PRIVATE calculus const fsutils grid mpimod timeutils)
target_compile_options(io PRIVATE ${FFLAGS})

if(hdf5)
  target_include_directories(io PUBLIC ${HDF5_INCLUDE_DIRS} ${HDF5_Fortran_INCLUDE_DIRS})
  target_compile_options(io PUBLIC ${HDF5_Fortran_DEFINITIONS})
  target_link_libraries(io PUBLIC h5fortran ${HDF5_Fortran_HL_LIBRARIES} ${HDF5_Fortran_LIBRARIES})
endif()
