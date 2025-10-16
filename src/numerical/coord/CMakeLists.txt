set_property(DIRECTORY PROPERTY  LABELS coord)

add_library(newton OBJECT newton.f90)
target_link_libraries(newton PRIVATE const)

add_library(spherical OBJECT spherical.f90)
target_link_libraries(spherical PRIVATE const)

add_library(dipole OBJECT dipole.f90)
target_link_libraries(dipole PRIVATE newton)

add_library(geomagnetic OBJECT geomagnetic.f90)
target_link_libraries(geomagnetic PRIVATE const)

