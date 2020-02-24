if(realbits EQUAL 32)
  message(STATUS "32-bit MUMPS interface")
  add_library(mumps_interface OBJECT mumps_real32.f90)
else()
  message(STATUS "64-bit MUMPS interface")
  add_library(mumps_interface OBJECT mumps_real64.f90)
endif()
target_include_directories(mumps_interface PUBLIC ${MUMPS_INCLUDE_DIRS})
target_link_libraries(mumps_interface PUBLIC ${MUMPS_LIBRARIES})

add_library(PDEelliptic OBJECT PDEelliptic.f90 elliptic2d.f90 elliptic3d.f90)
target_link_libraries(PDEelliptic PUBLIC const mumps_interface MPI::MPI_Fortran)

add_library(potential OBJECT
potential_comm_mumps.f90 potential_root.f90 potential_worker.f90
potential_mumps.f90 potential2d.f90
boundary_conditions/potentialBCs_mumps.f90)
target_link_libraries(potential PUBLIC calculus collision const grid mesh interp timeutils PDEelliptic mpimod reader MPI::MPI_Fortran)

if(NOT hdf5)
  return()
endif()

set(_potential_testfile ${CMAKE_CURRENT_BINARY_DIR}/test_potential2d.h5)
add_executable(test_potential2d test_potential2D.f90)
target_link_libraries(test_potential2d const PDEelliptic MPI::MPI_Fortran h5fortran)
add_test(NAME unit:potential2d COMMAND test_potential2d WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(unit:potential2d PROPERTIES TIMEOUT 30
                     FIXTURES_SETUP GemPot)

if(OctaveOK)
add_test(NAME unit:potential2d_octave
  COMMAND ${Octave_EXECUTABLE} -q --eval "test_potential2D('${_potential_testfile}')"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_tests_properties(unit:potential2d_octave PROPERTIES
  REQUIRED_FILES ${_potential_testfile}
  TIMEOUT 15
  FIXTURES_REQUIRED GemPot)
endif(OctaveOK)

if(PythonOK)
add_test(NAME unit:potential2d_python
  COMMAND ${Python3_EXECUTABLE}  ${CMAKE_CURRENT_SOURCE_DIR}/test_potential2d.py ${_potential_testfile})
set_tests_properties(unit:potential2d_python PROPERTIES
  REQUIRED_FILES ${_potential_testfile}
  TIMEOUT 10
  FIXTURES_REQUIRED GemPot)
endif(PythonOK)
