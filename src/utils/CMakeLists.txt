#--------magnetic field calculation executable---------------
configure_file(magcalc_cli.in.f90 magcalc_cli.f90 @ONLY)
add_executable(magcalc.bin)
target_sources(magcalc.bin PRIVATE magcalc.f90 ${CMAKE_CURRENT_BINARY_DIR}/magcalc_cli.f90
$<TARGET_OBJECTS:autogrid> $<TARGET_OBJECTS:grid_mpi>  $<TARGET_OBJECTS:grid>
$<TARGET_OBJECTS:mpimod>
$<TARGET_OBJECTS:meshobj> $<TARGET_OBJECTS:meshobj_cart> $<TARGET_OBJECTS:meshobj_dipole> $<TARGET_OBJECTS:spherical> $<TARGET_OBJECTS:newton> $<TARGET_OBJECTS:geomagnetic>
$<TARGET_OBJECTS:reader> $<TARGET_OBJECTS:io>
$<$<NOT:$<BOOL:${netcdf}>>:$<TARGET_OBJECTS:nc4fortran_dummy>>
)
target_link_libraries(magcalc.bin PRIVATE exe_frontend
sanity_check errors config timeutils const
hwloc_ifc hwloc_c
h5fortran::h5fortran
$<$<BOOL:${netcdf}>:nc4fortran::nc4fortran>
ffilesystem::filesystem
)
if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  set_target_properties(magcalc.bin PROPERTIES LINKER_LANGUAGE Fortran)
else()
  set_target_properties(magcalc.bin PROPERTIES LINKER_LANGUAGE CXX)
endif()
target_link_libraries(magcalc.bin PRIVATE MPI::MPI_Fortran
$<$<NOT:$<BOOL:${mpi}>>:MUMPS::MPISEQ>
)

# --- gemini3d.run Fortran front end
add_subdirectory(hwloc)

string(JOIN " " gemini_features
REALBITS:${realbits}
$<$<BOOL:${mpi}>:MPI>
$<$<BOOL:${glow}>:GLOW>
$<$<BOOL:${msis2}>:MSIS2>
$<$<BOOL:${hwm14}>:HWM14>
$<$<BOOL:${hdf5}>:HDF5>
$<$<BOOL:${netcdf}>:NETCDF4>
$<$<BOOL:${scotch}>:SCOTCH>
$<$<BOOL:${openmp}>:OPENMP>
$<$<BOOL:${HWLOC_FOUND}>:HWLOC>
)
configure_file(exe_frontend.in.f90 exe_frontend.f90.in @ONLY)
file(GENERATE OUTPUT exe_frontend.f90
INPUT ${CMAKE_CURRENT_BINARY_DIR}/exe_frontend.f90.in
)

add_library(exe_frontend OBJECT ${CMAKE_CURRENT_BINARY_DIR}/exe_frontend.f90)
target_link_libraries(exe_frontend PRIVATE config hwloc_ifc hwloc_c timeutils const
ffilesystem::filesystem
)

# --- gemini3d.run
add_executable(gemini3d.run gemini3d_run.f90
$<TARGET_OBJECTS:exe_frontend>
$<$<NOT:$<BOOL:${netcdf}>>:$<TARGET_OBJECTS:nc4fortran_dummy>>
)
target_link_libraries(gemini3d.run PRIVATE const config timeutils autogrid reader
h5fortran::h5fortran
$<$<BOOL:${netcdf}>:nc4fortran::nc4fortran>
hwloc_ifc hwloc_c
$<$<BOOL:${HWLOC_FOUND}>:HWLOC::HWLOC>
$<$<BOOL:${HDF5_HAVE_PARALLEL}>:MPI::MPI_Fortran>
ffilesystem::filesystem
)
if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  set_target_properties(gemini3d.run PROPERTIES LINKER_LANGUAGE Fortran)
else()
  set_target_properties(gemini3d.run PROPERTIES LINKER_LANGUAGE CXX)
endif()

# --- magcalc.run
add_executable(magcalc.run magcalc_run.f90
$<$<NOT:$<BOOL:${netcdf}>>:$<TARGET_OBJECTS:nc4fortran_dummy>>
)
target_link_libraries(magcalc.run PRIVATE exe_frontend autogrid reader config timeutils const hwloc_ifc hwloc_c
h5fortran::h5fortran
$<$<BOOL:${netcdf}>:nc4fortran::nc4fortran>
$<$<BOOL:${HWLOC_FOUND}>:HWLOC::HWLOC>
$<$<BOOL:${HDF5_HAVE_PARALLEL}>:MPI::MPI_Fortran>
ffilesystem::filesystem
)
if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  set_target_properties(magcalc.run PROPERTIES LINKER_LANGUAGE Fortran)
else()
  set_target_properties(magcalc.run PROPERTIES LINKER_LANGUAGE CXX)
endif()

# --- Gemini compare
add_library(assert OBJECT assert.f90)


set_property(TARGET gemini3d.run magcalc.bin magcalc.run
PROPERTY RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
)
set_property(TARGET gemini3d.run magcalc.bin magcalc.run
PROPERTY DEBUG_POSTFIX .debug
)
set_property(TARGET gemini3d.run magcalc.bin magcalc.run
PROPERTY RELWITHDEBINFO_POSTFIX .debug
)

# --- Install

install(TARGETS magcalc.bin gemini3d.run magcalc.run
EXPORT ${PROJECT_NAME}-targets
)
