set_property(DIRECTORY PROPERTY LABELS unit)

add_executable(test_compiler_vendor test_compiler_vendor.f90)
target_link_libraries(test_compiler_vendor PRIVATE gemini3d_sysinfo ffilesystem::filesystem)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  set(compiler_vendor_reg Intel)
else()
  set(compiler_vendor_reg ${CMAKE_Fortran_COMPILER_ID})
endif()

add_test(NAME CompilerVendor COMMAND test_compiler_vendor)
set_property(TEST CompilerVendor PROPERTY PASS_REGULAR_EXPRESSION ${compiler_vendor_reg})

add_executable(test_expand_envvar test_expand_envvar.f90)
target_link_libraries(test_expand_envvar PRIVATE gemini3d_config ffilesystem::filesystem)

add_test(NAME ExpandEnvvar COMMAND test_expand_envvar)
set_property(TEST ExpandEnvvar PROPERTY ENVIRONMENT test__gem=hello)

# apply compiler options
get_property(tgts DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)
foreach(t IN LISTS tgts)
  target_compile_options(${t} PRIVATE ${gemini3d_flags})
endforeach()
