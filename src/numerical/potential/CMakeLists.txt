# make mumps usage precision-agnostic
# this is needed because Mumps uses Fortran 95 syntax instead of Fortran 2003 syntax
configure_file(mumps_ifc.in.f90 mumps_ifc.f90 @ONLY)
add_library(mumps_interface OBJECT ${CMAKE_CURRENT_BINARY_DIR}/mumps_ifc.f90)
target_link_libraries(mumps_interface PRIVATE MUMPS::MUMPS)
# must be PUBLIC or gemini.bin[.exe] link will fail with GCC
# not OBJECT to avoid CMake 3.13 trouble

add_library(PDEelliptic OBJECT
PDEelliptic.f90 elliptic2d.f90 elliptic3d.f90)
target_link_libraries(PDEelliptic PRIVATE mumps_interface mpimod const MPI::MPI_Fortran)

add_library(potential OBJECT
potential_comm_mumps.f90 potential_root.f90 potential_worker.f90
potential_mumps.f90 potential2d.f90
boundary_conditions/potentialBCs_mumps.f90)
target_link_libraries(potential PRIVATE PDEelliptic calculus collision grid meshobj interp timeutils mpimod reader const MPI::MPI_Fortran)


set_targ_props(PDEelliptic potential)


if(hdf5 AND ${PROJECT_NAME}_BUILD_TESTING)

set(_potential_testfile ${CMAKE_CURRENT_BINARY_DIR}/test_potential2d.h5)

add_executable(test_potential2d test_potential2D.f90)
set_target_properties(test_potential2d PROPERTIES LABELS unit)
# MPI_C is needed here, else will not compile on MinGW with Ninja or Make.
# LINK_INTERFACE_MULTIPLICITY and IMPORTED_LINK_INTERFACE_MULTIPLICITY do not help here.
# adding INTERFACE to Mumps project didn't help either
# nor did adding them in mumps.cmake
# it seems to have to be right here.
if(mpi)
  target_link_libraries(test_potential2d PRIVATE PDEelliptic const h5fortran::h5fortran
    MUMPS::MUMPS SCALAPACK::SCALAPACK LAPACK::LAPACK MPI::MPI_Fortran)

  add_test(NAME unit:potential2
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${Ncpu} $<TARGET_FILE:test_potential2d> ${_potential_testfile})
else()
  target_link_libraries(test_potential2d PRIVATE PDEelliptic const h5fortran::h5fortran
    MUMPS::MUMPS LAPACK::LAPACK MPI::MPI_Fortran MUMPS::MPISEQ)
  # the reason for this is the MPI-1 Fortran 77 linking in MUMPS
  # this seemed to be an easy way to avoid symbol clashes

  add_test(NAME unit:potential2
    COMMAND $<TARGET_FILE:test_potential2d> ${_potential_testfile})

endif()
set_tests_properties(unit:potential2 PROPERTIES TIMEOUT 30
  FIXTURES_SETUP GemPot
  RESOURCE_LOCK cpu_mpi
  LABELS unit
  DEPENDS unit:HWLOC)


if(Python3_FOUND)

  add_test(NAME unit:potential2:python
    COMMAND ${Python3_EXECUTABLE}  ${CMAKE_CURRENT_SOURCE_DIR}/test_potential2d.py ${_potential_testfile})
  set_tests_properties(unit:potential2:python PROPERTIES
    REQUIRED_FILES ${_potential_testfile}
    TIMEOUT 10
    LABELS "unit;python"
    FIXTURES_REQUIRED GemPot)

endif(Python3_FOUND)

endif(hdf5 AND ${PROJECT_NAME}_BUILD_TESTING)
