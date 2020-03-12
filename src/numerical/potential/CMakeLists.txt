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
target_link_libraries(PDEelliptic PUBLIC const mpimod mumps_interface MPI::MPI_Fortran)

add_library(potential OBJECT
potential_comm_mumps.f90 potential_root.f90 potential_worker.f90
potential_mumps.f90 potential2d.f90
boundary_conditions/potentialBCs_mumps.f90)
target_link_libraries(potential PUBLIC calculus collision const grid mesh interp timeutils PDEelliptic mpimod reader MPI::MPI_Fortran)

if(NOT hdf5)
  return()
endif()

if(MINGW)
  message(VERBOSE "this test doesn't work with MS-MPI")
  return()
endif()

set(_potential_testfile ${CMAKE_CURRENT_BINARY_DIR}/test_potential2d.h5)
add_executable(test_potential2d test_potential2D.f90)
target_link_libraries(test_potential2d PRIVATE const PDEelliptic MPI::MPI_Fortran h5fortran)
add_test(NAME unit:potential2 COMMAND $<TARGET_FILE:test_potential2d>
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(unit:potential2 PROPERTIES TIMEOUT 30
  FIXTURES_SETUP GemPot)
win32_env(unit:potential2)

add_test(NAME unit:potential2:python
  COMMAND ${Python3_EXECUTABLE}  ${CMAKE_CURRENT_SOURCE_DIR}/test_potential2d.py ${_potential_testfile})
set_tests_properties(unit:potential2:python PROPERTIES
  REQUIRED_FILES ${_potential_testfile}
  TIMEOUT 10
  FIXTURES_REQUIRED GemPot
  DISABLED ${python_disabled})
