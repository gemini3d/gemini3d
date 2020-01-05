include(ExternalProject)

ExternalProject_Add(MUMPS
  GIT_REPOSITORY https://github.com/scivision/mumps.git
  GIT_TAG v5.2.1.4
  CMAKE_ARGS "-Darith=${arith}" "-Dparallel=true" "-Dmetis=${metis}" "-Dscotch=${scotch}" "-Dopenmp=false"
  INSTALL_COMMAND ""  # disables the install step for the external project
)

ExternalProject_Get_Property(MUMPS BINARY_DIR SOURCE_DIR)
set(MUMPS_BINARY_DIR ${BINARY_DIR})
set(MUMPS_SOURCE_DIR ${SOURCE_DIR})

# here we have to use a priori about MUMPS, since MUMPS won't build as ExernalProject at configure time,
# which is when find_package() is run
# An alternative (messy) is a make a "superproject" that uses both as external projects and calls CMake twice.
# This method below seems much preferable.
# Meson makes this much easier, and is a key reason to use Meson instead of CMake.

unset(MUMPS_LIBRARIES)
foreach(a ${arith})
  list(APPEND MUMPS_LIBRARIES ${a}mumps)
endforeach()
list(APPEND MUMPS_LIBRARIES mumps_common pord)

foreach(l ${MUMPS_LIBRARIES})
  add_library(${l} STATIC IMPORTED GLOBAL)
  add_dependencies(${l} MUMPS)
  set_target_properties(${l} PROPERTIES
    IMPORTED_LOCATION ${MUMPS_BINARY_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}${l}${CMAKE_STATIC_LIBRARY_SUFFIX})
endforeach()

set(MUMPS_INCLUDE_DIRS ${MUMPS_BINARY_DIR} ${MUMPS_SOURCE_DIR}/include)

foreach(a ${arith})
  target_link_libraries(${a}mumps INTERFACE mumps_common pord)
endforeach()
