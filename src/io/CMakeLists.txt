set_directory_properties(PROPERTIES LABELS io)

add_library(logging OBJECT logging.f90)

add_library(sanity_check OBJECT sanity_check.f90 $<TARGET_OBJECTS:errors>)
target_link_libraries(sanity_check PRIVATE const)

# --- error handling -- dump data to disk
add_library(errors OBJECT error_stop.f90)
target_link_libraries(errors PRIVATE const h5fortran::h5fortran)

# --- config
add_library(gemini3d_config OBJECT config.f90 config_nml.f90 config_ini.f90)
target_link_libraries(gemini3d_config PRIVATE const gemini3d_sysinfo
ffilesystem::filesystem
)

install(FILES ${PROJECT_BINARY_DIR}/include/gemini3d_config.mod TYPE INCLUDE)

# --- io
if(NOT netcdf)
  add_library(nc4fortran_dummy OBJECT $<$<NOT:$<BOOL:${netcdf}>>:${PROJECT_SOURCE_DIR}/src/vendor/nc4fortran_dummy.f90>)
  install(FILES ${PROJECT_BINARY_DIR}/include/nc4fortran.mod TYPE INCLUDE)
endif()

configure_file(output.in.f90 output.f90 @ONLY)
add_library(io OBJECT io.f90 ${CMAKE_CURRENT_BINARY_DIR}/output.f90
aurora.f90 aurora_raw.f90 aurora_hdf5.f90 aurora_ncdf.f90
cond.f90 cond_hdf5.f90
mag.f90 mag_hdf5.f90 mag_ncdf.f90 mag_raw.f90
plasma_input.f90 plasma_input_hdf5.f90 plasma_input_ncdf.f90 plasma_input_raw.f90
plasma_output.f90 plasma_output_hdf5.f90 plasma_output_ncdf.f90 plasma_output_raw.f90
milestone.f90
$<$<NOT:$<BOOL:${netcdf}>>:$<TARGET_OBJECTS:nc4fortran_dummy>>
)

target_link_libraries(io PRIVATE
calculus gemini3d_config const reader grid mpimod sanity_check timeutils
h5fortran::h5fortran
$<$<BOOL:${netcdf}>:nc4fortran::nc4fortran>
ffilesystem::filesystem
)

# --- reader
add_library(reader OBJECT reader.f90 reader_hdf5.f90 reader_ncdf.f90 reader_raw.f90
$<$<NOT:$<BOOL:${netcdf}>>:$<TARGET_OBJECTS:nc4fortran_dummy>>
)
target_link_libraries(reader PRIVATE const
h5fortran::h5fortran
$<$<BOOL:${netcdf}>:nc4fortran::nc4fortran>
ffilesystem::filesystem
)
