set_property(DIRECTORY PROPERTY  LABELS coord)

add_library(newton OBJECT newton.f90)
target_link_libraries(newton PRIVATE const)

add_library(spherical OBJECT spherical.f90)
target_link_libraries(spherical PRIVATE const)

add_library(geomagnetic OBJECT geomagnetic.f90)
target_link_libraries(geomagnetic PRIVATE const)

add_library(meshobj OBJECT meshobj.f90)
target_link_libraries(meshobj PRIVATE const geomagnetic spherical h5fortran::h5fortran)

add_library(meshobj_dipole OBJECT meshobj_dipole.f90 dipole_fns.f90
$<TARGET_OBJECTS:meshobj>
)
target_link_libraries(meshobj_dipole PRIVATE const newton spherical geomagnetic)

add_library(meshobj_cart OBJECT meshobj_cart.f90
$<TARGET_OBJECTS:meshobj>
)
target_link_libraries(meshobj_cart PRIVATE const spherical geomagnetic)
