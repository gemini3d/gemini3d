set_directory_properties(PROPERTIES LABELS "extra")


if(HAVE_BOOST_UT)

set(boost_ut_tests exepath hostname)

if(ffilesystem_locale)
  list(APPEND boost_ut_tests locale)
endif()

foreach(t IN LISTS boost_ut_tests)
  add_executable(test_${t}_cpp test_${t}.cpp)
  target_link_libraries(test_${t}_cpp PRIVATE ffilesystem)
  target_include_directories(test_${t}_cpp PRIVATE ${boost_ut_SOURCE_DIR}/include)
  target_compile_options(test_${t}_cpp PRIVATE ${ffilesystem_cpp_test_flags})

  add_test(NAME ${t}_cpp COMMAND test_${t}_cpp)
  set_tests_properties(${t}_cpp PROPERTIES LABELS "Cpp")
endforeach()

endif()

# --------- Fortran


if(ffilesystem_fortran)

set(extra_tests_fortran exepath filesystem_type)


foreach(t IN LISTS extra_tests_fortran)

  add_executable(test_${t}_fortran test_${t}.f90)
  target_link_libraries(test_${t}_fortran PRIVATE ffilesystem)
  target_compile_options(test_${t}_fortran PRIVATE ${ffilesystem_fortran_test_flags})
  set_property(TARGET test_${t}_fortran PROPERTY LINKER_LANGUAGE ${ffilesystem_linker_lang})

  add_test(NAME Fortran_${t} COMMAND test_${t}_fortran)
  set_tests_properties(Fortran_${t} PROPERTIES
  LABELS Fortran
  SKIP_RETURN_CODE 77
  WORKING_DIRECTORY $<TARGET_FILE_DIR:test_${t}_fortran>
  FIXTURES_SETUP ${t}
  )

endforeach()

target_compile_options(test_filesystem_type_fortran PRIVATE
$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wno-uninitialized>
)

set_tests_properties(Fortran_filesystem_type PROPERTIES
DISABLED $<AND:$<BOOL:${LINUX}>,$<NOT:$<BOOL:${ffilesystem_HAVE_LINUX_MAGIC}>>>
)

endif()
