set_directory_properties(PROPERTIES LABELS potential)

# make mumps usage precision-agnostic
# this is needed because Mumps uses Fortran 95 syntax instead of Fortran 2003 syntax
configure_file(mumps_ifc.in.f90 mumps_ifc.f90 @ONLY)
add_library(mumps_interface OBJECT ${CMAKE_CURRENT_BINARY_DIR}/mumps_ifc.f90)
target_link_libraries(mumps_interface PRIVATE MUMPS::MUMPS)

add_library(PDEelliptic OBJECT
PDEelliptic.f90 elliptic2d.f90 elliptic3d.f90)
target_link_libraries(PDEelliptic PRIVATE mpimod const mumps_interface MUMPS::MUMPS)

add_library(potential OBJECT
potential_comm_mumps.f90 potential_root.f90 potential_worker.f90
potential_mumps.f90 potential2d.f90
boundary_conditions/potentialBCs_mumps.f90)
target_link_libraries(potential PRIVATE PDEelliptic calculus collision grid meshobj efielddataobj interp timeutils mpimod reader const
inputdata efielddata
MPI::MPI_Fortran
)


if(NOT BUILD_TESTING)
  return()
endif()

set(_potential_testfile ${CMAKE_CURRENT_BINARY_DIR}/test_potential2d.h5)

add_executable(test_potential2d test_potential2D.f90 $<TARGET_OBJECTS:PDEelliptic>)
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

set(_cmd $<TARGET_FILE:test_potential2d> ${_potential_testfile})
if(mpi)
  list(PREPEND _cmd ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${Ncpu})
endif()
add_test(NAME unit:potential2 COMMAND ${_cmd})
set_tests_properties(unit:potential2 PROPERTIES TIMEOUT 30
FIXTURES_SETUP gemini_potential_fxt
RESOURCE_LOCK cpu_mpi
PROCESSORS ${Ncpu}
LABELS "unit"
DISABLED $<NOT:$<BOOL:${hdf5}>>
)

dll_test_path("h5fortran::h5fortran;MUMPS::MUMPS" unit:potential2)

if(python)

add_test(NAME unit:potential2:python
COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_potential2d.py ${_potential_testfile}
)
set_tests_properties(unit:potential2:python PROPERTIES
REQUIRED_FILES ${_potential_testfile}
TIMEOUT 10
LABELS "unit;python"
FIXTURES_REQUIRED gemini_potential_fxt
)

endif()
