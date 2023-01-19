set_property(DIRECTORY PROPERTY LABELS coord)

add_executable(newton_testdriver newton_testdriver.f90
$<TARGET_OBJECTS:meshobj>
)
target_link_libraries(newton_testdriver PRIVATE newton spherical geomagnetic meshobj_dipole const
h5fortran::h5fortran
)
add_test(NAME newton_testdriver COMMAND newton_testdriver)

add_executable(geomag2geog_testdriver geomag2geog_testdriver.f90)
target_link_libraries(geomag2geog_testdriver PRIVATE geomagnetic const)
add_test(NAME geomag2geog_testdriver COMMAND newton_testdriver)

add_executable(grid_testdriver grid_testdriver.f90)
target_link_libraries(grid_testdriver PRIVATE newton spherical geomagnetic meshobj_dipole meshobj const)
add_test(NAME grid_testdriver COMMAND grid_testdriver)

add_executable(fullgrid_dipole_testdriver fullgrid_dipole_testdriver.f90)
target_link_libraries(fullgrid_dipole_testdriver PRIVATE newton spherical geomagnetic meshobj_dipole meshobj const
ffilesystem::filesystem
)
add_test(NAME fullgrid_dipole_testdriver COMMAND fullgrid_dipole_testdriver)

add_executable(fullgrid_dipole_regen fullgrid_dipole_regen.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(fullgrid_dipole_regen PRIVATE reader newton spherical geomagnetic meshobj_dipole meshobj_cart meshobj mpimod autogrid const
h5fortran::h5fortran
ffilesystem::filesystem
)
add_test(NAME fullgrid_dipole_regen COMMAND fullgrid_dipole_regen)

add_executable(fullgrid_cartesian_testdriver fullgrid_cartesian_testdriver.f90)
target_link_libraries(fullgrid_cartesian_testdriver PRIVATE newton spherical geomagnetic meshobj_cart meshobj const
ffilesystem::filesystem
)
add_test(NAME fullgrid_cartesian_testdriver COMMAND fullgrid_cartesian_testdriver)

add_executable(fullgrid_dipole_testdriver_root fullgrid_dipole_testdriver_root.f90)
target_link_libraries(fullgrid_dipole_testdriver_root PRIVATE newton spherical geomagnetic meshobj_dipole meshobj const)
add_test(NAME fullgrid_dipole_testdriver_root COMMAND fullgrid_dipole_testdriver_root)

# --- test props
if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  set_property(TARGET fullgrid_dipole_regen fullgrid_dipole_testdriver fullgrid_cartesian_testdriver PROPERTY LINKER_LANGUAGE Fortran)
else()
  set_property(TARGET fullgrid_dipole_regen fullgrid_dipole_testdriver fullgrid_cartesian_testdriver PROPERTY LINKER_LANGUAGE CXX)
endif()

get_property(test_names DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY TESTS)

set_property(TEST ${test_names} PROPERTY LABELS unit)

if(host_ramGB LESS 8)
  set_property(TEST ${test_names} PROPERTY RESOURCE_LOCK cpu_ram)
endif()

dll_test_path("h5fortran::h5fortran;LAPACK::LAPACK;ffilesystem::filesystem" "${test_names}")
