set_property(DIRECTORY PROPERTY LABELS io)

add_library(logging OBJECT logging.f90)

add_library(sanity_check OBJECT sanity_check.f90 $<TARGET_OBJECTS:errors>)
target_link_libraries(sanity_check PRIVATE const)

# --- error handling -- dump data to disk
add_library(errors OBJECT error_stop.f90)
target_link_libraries(errors PRIVATE const h5fortran::h5fortran)

# --- config
add_library(gemini3d_config OBJECT config.f90 config_nml.f90 config_ini.f90)
target_link_libraries(gemini3d_config PRIVATE const timeutils gemini3d_sysinfo
ffilesystem::filesystem
)

install(FILES ${PROJECT_BINARY_DIR}/include/gemini3d_config.mod TYPE INCLUDE)

# --- io
add_library(atomic_file OBJECT atomic_file.f90 atomic_file.c)
add_library(restart_runtime OBJECT restart_runtime.f90)
target_link_libraries(restart_runtime PRIVATE atomic_file const gemini3d_config timeutils mpimod h5fortran::h5fortran MPI::MPI_Fortran)
configure_file(output.in.f90 output.f90 @ONLY)
add_library(io OBJECT io.f90 ${CMAKE_CURRENT_BINARY_DIR}/output.f90
aurora.f90 aurora_hdf5.f90
cond.f90 cond_hdf5.f90
mag.f90 mag_hdf5.f90
plasma_input.f90 plasma_input_hdf5.f90
plasma_output.f90 plasma_output_hdf5.f90
milestone.f90
$<TARGET_OBJECTS:grid>
)

target_link_libraries(io PRIVATE
calculus gemini3d_config const reader grid mpimod sanity_check timeutils gemini_work_def atomic_file restart_runtime
h5fortran::h5fortran
ffilesystem::filesystem
)

# --- io_nompi
add_library(io_nompi OBJECT io_nompi.f90)
target_link_libraries(io_nompi PRIVATE grid interp const timeutils reader
h5fortran::h5fortran
ffilesystem::filesystem
)

# --- reader
add_library(reader OBJECT reader.f90 reader_hdf5.f90)
target_link_libraries(reader PRIVATE const
h5fortran::h5fortran
ffilesystem::filesystem
)
