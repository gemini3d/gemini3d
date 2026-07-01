set_property(DIRECTORY PROPERTY LABELS Gemini3d:coord)

add_executable(newton_testdriver newton_testdriver.f90
$<TARGET_OBJECTS:meshobj>
)
target_link_libraries(newton_testdriver PRIVATE newton spherical geomagnetic meshobj_dipole dipole const
h5fortran::h5fortran
)
add_test(NAME gemini:NewtonTest COMMAND newton_testdriver)
hdf5_dll(gemini:NewtonTest)

add_executable(geomag2geog_testdriver geomag2geog_testdriver.f90)
target_link_libraries(geomag2geog_testdriver PRIVATE geomagnetic const)
add_test(NAME gemini:Geomag2Geog COMMAND geomag2geog_testdriver)
hdf5_dll(gemini:Geomag2Geog)

add_executable(grid_testdriver grid_testdriver.f90)
target_link_libraries(grid_testdriver PRIVATE newton spherical geomagnetic dipole meshobj_dipole meshobj const)
add_test(NAME gemini:GridTest COMMAND grid_testdriver)
hdf5_dll(gemini:GridTest)

add_executable(fullgrid_dipole_testdriver fullgrid_dipole_testdriver.f90)
target_link_libraries(fullgrid_dipole_testdriver PRIVATE newton spherical geomagnetic dipole meshobj_dipole meshobj const
ffilesystem::filesystem
)
add_test(NAME gemini:FullgridDipole COMMAND fullgrid_dipole_testdriver)
hdf5_dll(gemini:FullgridDipole)

add_executable(fullgrid_dipole_regen fullgrid_dipole_regen.f90
$<TARGET_OBJECTS:grid>
)
target_link_libraries(fullgrid_dipole_regen PRIVATE reader newton spherical geomagnetic dipole meshobj_dipole meshobj_cart meshobj mpimod autogrid const
h5fortran::h5fortran
ffilesystem::filesystem
)
add_test(NAME gemini:FullgridDipoleRegen COMMAND fullgrid_dipole_regen)
hdf5_dll(gemini:FullgridDipoleRegen)

add_executable(fullgrid_cartesian_testdriver fullgrid_cartesian_testdriver.f90)
target_link_libraries(fullgrid_cartesian_testdriver PRIVATE newton spherical geomagnetic meshobj_cart meshobj const
ffilesystem::filesystem
)
add_test(NAME gemini:FullgridCartesian COMMAND fullgrid_cartesian_testdriver)
hdf5_dll(gemini:FullgridCartesian)

add_executable(fullgrid_dipole_testdriver_root fullgrid_dipole_testdriver_root.f90)
target_link_libraries(fullgrid_dipole_testdriver_root PRIVATE newton spherical geomagnetic dipole meshobj_dipole meshobj const)
add_test(NAME gemini:FullgridDipoleRoot COMMAND fullgrid_dipole_testdriver_root)
hdf5_dll(gemini:FullgridDipoleRoot)

# --- test props
if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  set_property(TARGET fullgrid_dipole_regen fullgrid_dipole_testdriver fullgrid_cartesian_testdriver PROPERTY LINKER_LANGUAGE Fortran)
else()
  set_property(TARGET fullgrid_dipole_regen fullgrid_dipole_testdriver fullgrid_cartesian_testdriver PROPERTY LINKER_LANGUAGE CXX)
endif()

get_property(test_names DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY TESTS)

set_property(TEST ${test_names} PROPERTY LABELS Gemini3d:unit:coord)

if(host_ramGB LESS 8)
  set_property(TEST ${test_names} PROPERTY RESOURCE_LOCK cpu_ram)
endif()
