if(NOT BUILD_SHARED_LIBS)
  return()
endif()

if(CYGWIN)
  set(_p $<PATH:GET_FILENAME,ffilesystem>)
else()
  set(_p $<TARGET_FILE:ffilesystem>)
endif()

if(HAVE_BOOST_UT)

set(boost_ut_tests libpath)

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

add_executable(test_libpath_fortran test_libpath.f90)
target_link_libraries(test_libpath_fortran PRIVATE ffilesystem)
target_compile_options(test_libpath_fortran PRIVATE ${ffilesystem_fortran_test_flags})
set_property(TARGET test_libpath_fortran PROPERTY LINKER_LANGUAGE ${ffilesystem_linker_lang})

add_test(NAME Fortran_libpath COMMAND test_libpath_fortran)

set_tests_properties(Fortran_libpath PROPERTIES
PASS_REGULAR_EXPRESSION "${_p}"
LABELS "Fortran"
${_dis}
)

endif()
