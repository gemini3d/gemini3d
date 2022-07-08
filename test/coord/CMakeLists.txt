set_directory_properties(PROPERTIES LABELS coord)

add_executable(newton_testdriver newton_testdriver.f90)
target_link_libraries(newton_testdriver PRIVATE newton spherical geomagnetic meshobj_dipole meshobj const
$<$<BOOL:${HDF5_HAVE_PARALLEL}>:MPI::MPI_Fortran>
)
add_test(NAME newton_testdriver COMMAND newton_testdriver)

add_executable(geomag2geog_testdriver geomag2geog_testdriver.f90)
target_link_libraries(geomag2geog_testdriver PRIVATE geomagnetic const)
add_test(NAME geomag2geog_testdriver COMMAND newton_testdriver)

add_executable(grid_testdriver grid_testdriver.f90)
target_link_libraries(grid_testdriver PRIVATE newton spherical geomagnetic meshobj_dipole meshobj const
$<$<BOOL:${HDF5_HAVE_PARALLEL}>:MPI::MPI_Fortran>
)
add_test(NAME grid_testdriver COMMAND grid_testdriver)

add_executable(fullgrid_dipole_testdriver fullgrid_dipole_testdriver.f90)
target_link_libraries(fullgrid_dipole_testdriver PRIVATE newton spherical geomagnetic meshobj_dipole meshobj const
$<$<BOOL:${HDF5_HAVE_PARALLEL}>:MPI::MPI_Fortran>
ffilesystem::filesystem
)
add_test(NAME fullgrid_dipole_testdriver COMMAND fullgrid_dipole_testdriver)

add_executable(fullgrid_dipole_regen fullgrid_dipole_regen.f90
)
target_link_libraries(fullgrid_dipole_regen PRIVATE grid reader newton spherical geomagnetic meshobj_dipole meshobj const
h5fortran::h5fortran
$<$<BOOL:${HDF5_HAVE_PARALLEL}>:MPI::MPI_Fortran>
ffilesystem::filesystem
)
add_test(NAME fullgrid_dipole_regen COMMAND fullgrid_dipole_regen)

add_executable(fullgrid_cartesian_testdriver fullgrid_cartesian_testdriver.f90)
target_link_libraries(fullgrid_cartesian_testdriver PRIVATE newton spherical geomagnetic meshobj_cart meshobj const
$<$<BOOL:${HDF5_HAVE_PARALLEL}>:MPI::MPI_Fortran>
ffilesystem::filesystem
)
add_test(NAME fullgrid_cartesian_testdriver COMMAND fullgrid_cartesian_testdriver)

add_executable(fullgrid_dipole_testdriver_root fullgrid_dipole_testdriver_root.f90)
target_link_libraries(fullgrid_dipole_testdriver_root PRIVATE newton spherical geomagnetic meshobj_dipole meshobj const
$<$<BOOL:${HDF5_HAVE_PARALLEL}>:MPI::MPI_Fortran>
)
add_test(NAME fullgrid_dipole_testdriver_root COMMAND fullgrid_dipole_testdriver_root)

# --- test props
if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  set_target_properties(fullgrid_dipole_regen fullgrid_dipole_testdriver fullgrid_cartesian_testdriver
  PROPERTIES LINKER_LANGUAGE Fortran
  )
else()
  set_target_properties(fullgrid_dipole_regen fullgrid_dipole_testdriver fullgrid_cartesian_testdriver
  PROPERTIES LINKER_LANGUAGE CXX
  )
endif()

get_property(test_names DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY TESTS)

set_tests_properties(${test_names} PROPERTIES
TIMEOUT 90
LABELS unit
)

if(host_ramGB LESS 8)
  set_tests_properties(${test_names} PROPERTIES RESOURCE_LOCK cpu_ram)
endif()

dll_test_path("h5fortran::h5fortran;LAPACK::LAPACK" "${test_names}")
