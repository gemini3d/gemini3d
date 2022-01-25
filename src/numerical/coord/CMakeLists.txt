set_directory_properties(PROPERTIES LABELS coord)

add_library(newton OBJECT newton.f90)
target_link_libraries(newton PRIVATE const)

add_library(spherical OBJECT spherical.f90)
target_link_libraries(spherical PRIVATE const)

add_library(geomagnetic OBJECT geomagnetic.f90)
target_link_libraries(geomagnetic PRIVATE const)

add_library(meshobj OBJECT meshobj.f90)
target_link_libraries(meshobj PRIVATE const geomagnetic spherical h5fortran::h5fortran)

add_library(meshobj_dipole OBJECT meshobj_dipole.f90 dipole_fns.f90)
target_link_libraries(meshobj_dipole PRIVATE const meshobj newton spherical geomagnetic)

add_library(meshobj_cart OBJECT meshobj_cart.f90)
target_link_libraries(meshobj_cart PRIVATE const meshobj spherical geomagnetic)

# --- coord self-tests

if(NOT BUILD_TESTING)
  return()
endif()

add_executable(newton_testdriver newton_testdriver.f90)
target_link_libraries(newton_testdriver PRIVATE newton spherical geomagnetic meshobj_dipole meshobj const
$<$<BOOL:${HDF5_HAVE_PARALLEL}>:MPI::MPI_Fortran>
)

add_executable(geomag2geog_testdriver geomag2geog_testdriver.f90)
target_link_libraries(geomag2geog_testdriver PRIVATE geomagnetic const)

add_executable(grid_testdriver grid_testdriver.f90)
target_link_libraries(grid_testdriver PRIVATE newton spherical geomagnetic meshobj_dipole meshobj const
$<$<BOOL:${HDF5_HAVE_PARALLEL}>:MPI::MPI_Fortran>
)

add_executable(fullgrid_dipole_testdriver fullgrid_dipole_testdriver.f90)
target_link_libraries(fullgrid_dipole_testdriver PRIVATE newton spherical geomagnetic meshobj_dipole meshobj const
$<$<BOOL:${HDF5_HAVE_PARALLEL}>:MPI::MPI_Fortran>
ffilesystem::filesystem ${lib_filesystem}
)
if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  set_target_properties(fullgrid_dipole_testdriver PROPERTIES LINKER_LANGUAGE Fortran)
else()
  set_target_properties(fullgrid_dipole_testdriver PROPERTIES LINKER_LANGUAGE CXX)
endif()

add_executable(fullgrid_cartesian_testdriver fullgrid_cartesian_testdriver.f90)
target_link_libraries(fullgrid_cartesian_testdriver PRIVATE newton spherical geomagnetic meshobj_cart meshobj const
$<$<BOOL:${HDF5_HAVE_PARALLEL}>:MPI::MPI_Fortran>
ffilesystem::filesystem ${lib_filesystem}
)
if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  set_target_properties(fullgrid_cartesian_testdriver PROPERTIES LINKER_LANGUAGE Fortran)
else()
  set_target_properties(fullgrid_cartesian_testdriver PROPERTIES LINKER_LANGUAGE CXX)
endif()

add_executable(fullgrid_dipole_testdriver_root fullgrid_dipole_testdriver_root.f90)
target_link_libraries(fullgrid_dipole_testdriver_root PRIVATE newton spherical geomagnetic meshobj_dipole meshobj const
$<$<BOOL:${HDF5_HAVE_PARALLEL}>:MPI::MPI_Fortran>
)

foreach(t newton_testdriver geomag2geog_testdriver grid_testdriver fullgrid_dipole_testdriver fullgrid_cartesian_testdriver fullgrid_dipole_testdriver_root)

  add_test(NAME "unit:${t}" COMMAND $<TARGET_FILE:${t}>)

  set_tests_properties(unit:${t} PROPERTIES
  TIMEOUT 30
  LABELS "unit"
  )

endforeach()

set_tests_properties(unit:fullgrid_dipole_testdriver unit:fullgrid_dipole_testdriver_root unit:fullgrid_cartesian_testdriver unit:grid_testdriver unit:newton_testdriver PROPERTIES
ENVIRONMENT $<$<BOOL:${test_dll_path}>:"PATH=${test_dll_path}">
)


if(host_ramGB LESS 8)
  set_tests_properties(unit:fullgrid_cartesian_testdriver unit:fullgrid_dipole_testdriver unit:fullgrid_dipole_testdriver_root PROPERTIES
  RESOURCE_LOCK cpu_ram
  )
endif()
