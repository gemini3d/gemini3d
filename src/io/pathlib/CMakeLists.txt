set_directory_properties(PROPERTIES LABELS pathlib)

add_library(pathlib OBJECT pathlib.f90)

if(WIN32)
  target_sources(pathlib PRIVATE pathlib_windows.f90)
else()
  target_sources(pathlib PRIVATE pathlib_unix.f90)
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  target_sources(pathlib PRIVATE pathlib_intel.f90)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  target_sources(pathlib PRIVATE pathlib_gcc.f90)
else()
  message(STATUS "pathlib_dummy: better to write pathlib_gcc-like routine for compiler ${CMAKE_Fortran_COMPILER_ID}")
  target_sources(pathlib PRIVATE pathlib_dummy.f90)
endif()

if(BUILD_TESTING)
  add_executable(test_pathlib test_pathlib.f90)
  target_link_libraries(test_pathlib PRIVATE pathlib)

  add_test(NAME unit:pathlib COMMAND $<TARGET_FILE:test_pathlib>)
  set_tests_properties(unit:pathlib PROPERTIES
    LABELS unit
    TIMEOUT 10)
endif()
