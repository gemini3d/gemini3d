add_library(gemini3d_sysinfo OBJECT sys.f90)
target_link_libraries(gemini3d_sysinfo PRIVATE ffilesystem::filesystem)

#--------magnetic field calculation executable---------------
configure_file(magcalc_cli.in.f90 magcalc_cli.f90 @ONLY)
add_executable(magcalc.bin
magcalc.f90 ${CMAKE_CURRENT_BINARY_DIR}/magcalc_cli.f90
$<TARGET_OBJECTS:autogrid> $<TARGET_OBJECTS:grid>
$<TARGET_OBJECTS:mpimod>
$<TARGET_OBJECTS:meshobj> $<TARGET_OBJECTS:meshobj_cart> $<TARGET_OBJECTS:meshobj_dipole> $<TARGET_OBJECTS:spherical> $<TARGET_OBJECTS:newton> $<TARGET_OBJECTS:geomagnetic>
$<TARGET_OBJECTS:reader> $<TARGET_OBJECTS:io> $<TARGET_OBJECTS:interp>
)
target_link_libraries(magcalc.bin PRIVATE exe_frontend
sanity_check errors gemini3d_config timeutils const gemini3d_sysinfo
h5fortran::h5fortran
ffilesystem::filesystem
)
# avoid double-linking out of order, link failure
target_link_libraries(magcalc.bin PRIVATE $<$<NOT:$<BOOL:${HDF5_HAVE_PARALLEL}>>:MPI::MPI_Fortran>)

# --- gemini3d.run Fortran front end

string(JOIN " " gemini_features
REALBITS:${realbits}
MPI
MSIS2
)
if(glow)
  string(APPEND gemini_features " GLOW")
endif()
if(hwm14)
  string(APPEND gemini_features " HWM14")
endif()
configure_file(exe_frontend.in.f90 exe_frontend.f90 @ONLY)

add_library(exe_frontend OBJECT cpu_count.cpp ${CMAKE_CURRENT_BINARY_DIR}/exe_frontend.f90)
target_compile_features(exe_frontend PRIVATE cxx_std_11)
target_compile_definitions(exe_frontend PRIVATE $<$<BOOL:${MSVC}>:_CRT_SECURE_NO_WARNINGS>)
target_link_libraries(exe_frontend PRIVATE gemini3d_config timeutils const
ffilesystem::filesystem
)

# --- gemini3d.run
add_executable(gemini3d.run gemini3d_run.f90
$<TARGET_OBJECTS:exe_frontend>
)
target_link_libraries(gemini3d.run PRIVATE
const gemini3d_config timeutils autogrid reader gemini3d_sysinfo
h5fortran::h5fortran
ffilesystem::filesystem
)

# --- magcalc.run
add_executable(magcalc.run magcalc_run.f90
$<TARGET_OBJECTS:exe_frontend>
)
target_link_libraries(magcalc.run PRIVATE
autogrid reader gemini3d_config timeutils const gemini3d_sysinfo
h5fortran::h5fortran
ffilesystem::filesystem
)


if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  set_property(TARGET gemini3d.run magcalc.run magcalc.bin PROPERTY LINKER_LANGUAGE Fortran)
else()
  set_property(TARGET gemini3d.run magcalc.run magcalc.bin PROPERTY LINKER_LANGUAGE CXX)
endif()

# --- Gemini compare
add_library(assert OBJECT assert.f90)


set_property(TARGET gemini3d.run magcalc.bin magcalc.run PROPERTY RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
set_property(TARGET gemini3d.run magcalc.bin magcalc.run PROPERTY DEBUG_POSTFIX .debug)
set_property(TARGET gemini3d.run magcalc.bin magcalc.run PROPERTY RELWITHDEBINFO_POSTFIX .debug)

# --- Install

install(TARGETS magcalc.bin gemini3d.run magcalc.run
EXPORT ${PROJECT_NAME}-targets
)
