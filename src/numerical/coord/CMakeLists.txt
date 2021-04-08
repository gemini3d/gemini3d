add_library(newton OBJECT newton.f90)
target_link_libraries(newton PRIVATE const)

add_library(meshobj OBJECT meshobj.f90)
target_link_libraries(meshobj PRIVATE const)

add_library(meshobj_dipole OBJECT meshobj_dipole.f90 dipole_fns.f90)
target_link_libraries(meshobj_dipole PRIVATE const meshobj newton)

add_executable(newton_testdriver newton_testdriver.f90)
target_link_libraries(newton_testdriver PRIVATE const newton meshobj meshobj_dipole)

add_executable(geomag2geog_testdriver geomag2geog_testdriver.f90)
target_link_libraries(geomag2geog_testdriver PRIVATE const newton meshobj meshobj_dipole)

add_executable(grid_testdriver grid_testdriver.f90)
target_link_libraries(grid_testdriver PRIVATE const newton meshobj meshobj_dipole)

add_executable(fullgrid_testdriver fullgrid_testdriver.f90)
target_link_libraries(fullgrid_testdriver PRIVATE const newton meshobj meshobj_dipole)
