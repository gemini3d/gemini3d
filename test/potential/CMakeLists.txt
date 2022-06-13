set_directory_properties(PROPERTIES LABELS potential)

# --- potential 2D test

set(_potential2d_testfile ${CMAKE_CURRENT_BINARY_DIR}/test_potential2d.h5)

add_executable(test_potential2d test_potential2d.f90 $<TARGET_OBJECTS:PDEelliptic>)
set_target_properties(test_potential2d PROPERTIES LABELS "unit")
target_link_libraries(test_potential2d PRIVATE const h5fortran::h5fortran MUMPS::MUMPS
$<$<NOT:$<BOOL:${mpi}>>:MUMPS::MPISEQ>
"$<$<BOOL:${MUMPS_HAVE_Scotch}>:Scotch::Scotch;METIS::METIS>"
)
# MPISEQ is the MPI-1 Fortran 77 linking in MUMPS--avoids symbol clashes

# LINK_INTERFACE_MULTIPLICITY and IMPORTED_LINK_INTERFACE_MULTIPLICITY do not help here.
# adding INTERFACE to Mumps project didn't help either
# nor did adding them in mumps.cmake
# it seems to have to be right here.
# test_potential2d has been an intermittent source of link-order issues,
# where mpi gets linked in the middle instead of at the end.
# the TARGET_OBJECTS etc in this directory are a careful effort to avoid this problem

set(_cmd $<TARGET_FILE:test_potential2d> ${_potential2d_testfile})
if(mpi)
  list(PREPEND _cmd ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${Ncpu})
endif()
add_test(NAME potential2d COMMAND ${_cmd})
set_tests_properties(potential2d PROPERTIES TIMEOUT 30
FIXTURES_SETUP gemini_potential_fxt
RESOURCE_LOCK cpu_mpi
PROCESSORS ${Ncpu}
LABELS "unit"
DISABLED $<NOT:$<BOOL:${hdf5}>>
)

dll_test_path("h5fortran::h5fortran;MUMPS::MUMPS" potential2d)

if(H5PY_FOUND)

add_test(NAME potential2d:python
COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_potential2d.py ${_potential2d_testfile}
)
set_tests_properties(potential2d:python PROPERTIES
REQUIRED_FILES ${_potential2d_testfile}
TIMEOUT 10
LABELS "unit;python"
FIXTURES_REQUIRED gemini_potential_fxt
)

endif()


# --- potential 3D test

set(_potential3d_testfile ${CMAKE_CURRENT_BINARY_DIR}/test_potential3d.h5)

add_executable(test_potential3d test_potential3d.f90 $<TARGET_OBJECTS:PDEelliptic>)
set_target_properties(test_potential3d PROPERTIES LABELS "unit")
target_link_libraries(test_potential3d PRIVATE const h5fortran::h5fortran MUMPS::MUMPS
$<$<NOT:$<BOOL:${mpi}>>:MUMPS::MPISEQ>
"$<$<BOOL:${MUMPS_HAVE_Scotch}>:Scotch::Scotch;METIS::METIS>"
)

set(_cmd $<TARGET_FILE:test_potential3d> ${_potential3d_testfile})
if(mpi)
  list(PREPEND _cmd ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${Ncpu})
endif()
add_test(NAME potential3d COMMAND ${_cmd})
set_tests_properties(potential3d PROPERTIES TIMEOUT 30
FIXTURES_SETUP gemini_potential_fxt
RESOURCE_LOCK cpu_mpi
PROCESSORS ${Ncpu}
LABELS "unit"
DISABLED $<NOT:$<BOOL:${hdf5}>>
)

dll_test_path("h5fortran::h5fortran;MUMPS::MUMPS" potential3d)
