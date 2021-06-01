set_directory_properties(PROPERTIES LABELS coord)

add_library(newton OBJECT newton.f90)
target_link_libraries(newton PRIVATE const)

add_library(spherical OBJECT spherical.f90)
target_link_libraries(spherical PRIVATE const)

add_library(geomagnetic OBJECT geomagnetic.f90)
target_link_libraries(geomagnetic PRIVATE const)

add_library(meshobj OBJECT meshobj.f90)
target_link_libraries(meshobj PRIVATE const geomagnetic h5fortran::h5fortran)

add_library(meshobj_dipole OBJECT meshobj_dipole.f90 dipole_fns.f90)
target_link_libraries(meshobj_dipole PRIVATE const meshobj newton spherical geomagnetic)

add_library(meshobj_cart OBJECT meshobj_cart.f90)
target_link_libraries(meshobj_cart PRIVATE const meshobj spherical geomagnetic)

set_targ_props(meshobj meshobj_cart meshobj_dipole newton spherical geomagnetic)

# --- coord self-tests

if(${PROJECT_NAME}_BUILD_TESTING)
add_executable(newton_testdriver newton_testdriver.f90)
target_link_libraries(newton_testdriver PRIVATE const newton meshobj meshobj_dipole newton spherical geomagnetic)

add_executable(geomag2geog_testdriver geomag2geog_testdriver.f90)
target_link_libraries(geomag2geog_testdriver PRIVATE const geomagnetic)

add_executable(grid_testdriver grid_testdriver.f90)
target_link_libraries(grid_testdriver PRIVATE const meshobj meshobj_dipole newton spherical geomagnetic)

add_executable(fullgrid_dipole_testdriver fullgrid_dipole_testdriver.f90)
target_link_libraries(fullgrid_dipole_testdriver PRIVATE const meshobj meshobj_dipole newton spherical geomagnetic)

add_executable(fullgrid_cartesian_testdriver fullgrid_cartesian_testdriver.f90)
target_link_libraries(fullgrid_cartesian_testdriver PRIVATE const meshobj meshobj_cart newton spherical geomagnetic)

add_executable(fullgrid_dipole_testdriver_root fullgrid_dipole_testdriver_root.f90)
target_link_libraries(fullgrid_dipole_testdriver_root PRIVATE const meshobj meshobj_dipole newton spherical geomagnetic)

foreach(t newton_testdriver geomag2geog_testdriver grid_testdriver fullgrid_dipole_testdriver fullgrid_cartesian_testdriver fullgrid_dipole_testdriver_root)

  add_test(NAME ${t} COMMAND $<TARGET_FILE:${t}>)

  set_tests_properties(${t} PROPERTIES
    TIMEOUT 10
    LABELS "unit"
    )

endforeach()

endif()
