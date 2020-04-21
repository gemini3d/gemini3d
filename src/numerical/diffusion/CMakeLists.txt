add_library(PDEparabolic OBJECT PDEparabolic.f90)
target_link_libraries(PDEparabolic PRIVATE const gbsv ${LAPACK_LIBRARIES})
target_include_directories(PDEparabolic PRIVATE ${LAPACK_INCLUDE_DIRS})

add_library(diffusion OBJECT diffusion.f90)
target_link_libraries(diffusion PRIVATE const config grid mesh PDEparabolic)

if(NOT hdf5)
  return()
endif()

# --- tests

set(_diffstem ${CMAKE_CURRENT_BINARY_DIR}/test_diffusion1d)

add_executable(test_diffusion1d test_diffusion1D.f90)
# lapack needed for MacOS
target_link_libraries(test_diffusion1d PRIVATE PDEparabolic gbsv const
  ${LAPACK_LIBRARIES} h5fortran::h5fortran)



add_test(NAME unit:diffusion1:dirichlet COMMAND $<TARGET_FILE:test_diffusion1d> 1 "${_diffstem}_dirichlet.h5")
set_tests_properties(unit:diffusion1:dirichlet PROPERTIES
  TIMEOUT 5
  FIXTURES_SETUP GemDiffDirich)
win32_env(unit:diffusion1:dirichlet)

add_test(NAME unit:diffusion1:neumann COMMAND $<TARGET_FILE:test_diffusion1d> 0 "${_diffstem}_neumann.h5")
set_tests_properties(unit:diffusion1:neumann PROPERTIES
  TIMEOUT 5
  FIXTURES_SETUP GemDiffNeumann)
win32_env(unit:diffusion1:neumann)

if(python_ok)
  add_test(NAME unit:diffusion1:dirichlet:python
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_diffusion1D.py "${_diffstem}_dirichlet.h5")
  set_tests_properties(unit:diffusion1:dirichlet:python PROPERTIES
    TIMEOUT 10
    REQUIRED_FILES "${_diffstem}_dirichlet.h5"
    FIXTURES_REQUIRED GemDiffDirich)

  add_test(NAME unit:diffusion1:neumann:python
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_diffusion1D.py "${_diffstem}_neumann.h5")
  set_tests_properties(unit:diffusion1:neumann:python PROPERTIES
    TIMEOUT 10
    REQUIRED_FILES "${_diffstem}_neumann.h5"
    FIXTURES_REQUIRED GemDiffNeumann)
endif(python_ok)