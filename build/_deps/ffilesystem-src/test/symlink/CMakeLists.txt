
if(HAVE_BOOST_UT)

set(boost_ut_tests is_symlink symlink)

foreach(t IN LISTS boost_ut_tests)
  add_executable(test_${t}_cpp test_${t}.cpp)
  target_link_libraries(test_${t}_cpp PRIVATE ffilesystem)
  target_include_directories(test_${t}_cpp PRIVATE ${boost_ut_SOURCE_DIR}/include)
  target_compile_options(test_${t}_cpp PRIVATE ${ffilesystem_cpp_test_flags})

  add_test(NAME ${t}_cpp COMMAND test_${t}_cpp)
  set_tests_properties(${t}_cpp PROPERTIES LABELS "Cpp")
endforeach()

endif()


if(ffilesystem_fortran)

foreach(t IN ITEMS is_symlink symlink)
  add_executable(test_${t}_fortran test_${t}.f90)
  target_link_libraries(test_${t}_fortran PRIVATE ffilesystem)
  target_compile_options(test_${t}_fortran PRIVATE ${ffilesystem_fortran_test_flags})
  set_property(TARGET test_${t}_fortran PROPERTY LINKER_LANGUAGE ${ffilesystem_linker_lang})
endforeach()

add_test(NAME Fortran_is_symlink COMMAND test_is_symlink_fortran)
set_tests_properties(Fortran_is_symlink PROPERTIES
FIXTURES_SETUP is_symlink
DISABLED $<BOOL:${sym_ret}>
)

add_test(NAME Fortran_symlink COMMAND test_symlink_fortran)
set_tests_properties(Fortran_symlink PROPERTIES
FIXTURES_SETUP symlink
FIXTURES_REQUIRED is_symlink
SKIP_RETURN_CODE 77
)

endif()
