if(HAVE_CXX_FILESYSTEM)
  set(cli_regex "Backend: <filesystem>")
else()
  set(cli_regex "Backend: C")
endif()


if(TARGET filesystem_cli)

add_test(NAME Fortran_CLI
COMMAND ${CMAKE_COMMAND} -Dexe:FILEPATH=$<TARGET_FILE:filesystem_cli> -P ${CMAKE_CURRENT_SOURCE_DIR}/stdin_nul.cmake
)
set_tests_properties(Fortran_CLI PROPERTIES
LABELS "Fortran"
PASS_REGULAR_EXPRESSION ${cli_regex}
)

endif()


if(TARGET fs_cli)

add_test(NAME Cpp_CLI
COMMAND ${CMAKE_COMMAND} -Dexe:FILEPATH=$<TARGET_FILE:fs_cli> -P ${CMAKE_CURRENT_SOURCE_DIR}/stdin_nul.cmake
)

set_tests_properties(Cpp_CLI PROPERTIES
LABELS "Cpp"
PASS_REGULAR_EXPRESSION ${cli_regex}
)

add_test(NAME CppCLInoLeak
COMMAND ${CMAKE_COMMAND} -Dexe:FILEPATH=$<TARGET_FILE:fs_cli> -P ${CMAKE_CURRENT_SOURCE_DIR}/stdin.cmake
)

set_tests_properties(CppCLInoLeak PROPERTIES
DISABLED $<AND:$<BOOL:${LINUX}>,$<NOT:$<BOOL:${ffilesystem_HAVE_LINUX_MAGIC}>>>
LABELS "Cpp"
)

endif()
