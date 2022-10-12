set_property(DIRECTORY PROPERTY LABELS utils)

add_executable(test_compiler_vendor test_compiler_vendor.f90)
target_link_libraries(test_compiler_vendor PRIVATE
gemini3d_sysinfo
ffilesystem::filesystem
)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  set(compiler_vendor_reg Intel)
else()
  set(compiler_vendor_reg ${CMAKE_Fortran_COMPILER_ID})
endif()

add_test(NAME CompilerVendor COMMAND test_compiler_vendor)
set_tests_properties(CompilerVendor PROPERTIES
PASS_REGULAR_EXPRESSION ${compiler_vendor_reg}
LABELS unit
)

add_executable(test_expand_envvar test_expand_envvar.f90)
target_link_libraries(test_expand_envvar PRIVATE gemini3d_config
ffilesystem::filesystem
)

add_test(NAME ExpandEnvvar COMMAND test_expand_envvar)
set_tests_properties(ExpandEnvvar PROPERTIES
ENVIRONMENT test__gem=hello
LABELS unit
)

# --- test props

get_property(test_names DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY TESTS)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  set_target_properties(${test_names} PROPERTIES LINKER_LANGUAGE Fortran)
else()
  set_target_properties(${test_names} PROPERTIES LINKER_LANGUAGE CXX)
endif()

dll_test_path(ffilesystem::filesystem "${test_names}")
