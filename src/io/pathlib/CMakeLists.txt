set_directory_properties(PROPERTIES LABELS pathlib)

add_library(pathlib OBJECT pathlib.f90
$<IF:$<BOOL:${WIN32}>,pathlib_windows.f90,pathlib_unix.f90>
$<$<Fortran_COMPILER_ID:Intel,IntelLLVM>:pathlib_intel.f90>
$<$<Fortran_COMPILER_ID:GNU>:pathlib_gcc.f90>
$<$<NOT:$<Fortran_COMPILER_ID:GNU,Intel,IntelLLVM>>:pathlib_dummy.f90>
)

if(BUILD_TESTING)
  add_executable(test_pathlib test_pathlib.f90)
  target_link_libraries(test_pathlib PRIVATE pathlib)

  add_test(NAME unit:pathlib COMMAND $<TARGET_FILE:test_pathlib>)
  set_tests_properties(unit:pathlib PROPERTIES
    LABELS unit
    TIMEOUT 10)
endif()
