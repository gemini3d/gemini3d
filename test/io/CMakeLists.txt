set_property(DIRECTORY PROPERTY LABELS io)

set(_test_file ${PROJECT_SOURCE_DIR}/test/config/config_example.nml)

add_executable(test_namelist test_nml.f90)
add_test(NAME Fortran_namelist COMMAND test_namelist ${_test_file})
set_tests_properties(Fortran_namelist PROPERTIES
REQUIRED_FILES ${_test_file}
LABELS unit
)


if(f2018assumed_rank)
  add_library(debug_dump OBJECT debug_dump.f90)
  target_link_libraries(debug_dump PRIVATE const h5fortran::h5fortran)
  set_target_properties(debug_dump PROPERTIES EXCLUDE_FROM_ALL true)
endif()

add_executable(test_dump_write test_dump_write.f90)
target_link_libraries(test_dump_write PRIVATE errors const
h5fortran::h5fortran
)

add_executable(test_dump_read test_dump_read.f90)
target_link_libraries(test_dump_read PRIVATE h5fortran::h5fortran)

foreach(t root worker input)
  set(f ${CMAKE_CURRENT_BINARY_DIR}/dump_${t}.h5)

  add_test(NAME gemini:dump:${t} COMMAND test_dump_write ${t} ${f})
  set_tests_properties(gemini:dump:${t} PROPERTIES
  FIXTURES_SETUP gemini:dump_fxt
  WILL_FAIL true
  LABELS unit
  )

  add_test(NAME gemini:dump:check:${t} COMMAND test_dump_read testing:${t} ${f})
  set_tests_properties(gemini:dump:check:${t} PROPERTIES
  FIXTURES_REQUIRED gemini:dump_fxt
  REQUIRED_FILES ${f}
  LABELS unit
  )

  dll_test_path("h5fortran::h5fortran;HDF5::HDF5" "gemini:dump:${t};gemini:dump:check:${t}")
endforeach()
