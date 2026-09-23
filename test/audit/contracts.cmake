add_executable(audit_partition partition_contract.f90)
target_link_libraries(audit_partition PRIVATE autogrid)
add_test(NAME audit:partition_contract COMMAND audit_partition)
add_test(NAME audit:partition_cmake COMMAND ${CMAKE_COMMAND}
  -Dsource=${PROJECT_SOURCE_DIR} -P ${CMAKE_CURRENT_LIST_DIR}/partition_contract.cmake)
set_tests_properties(audit:partition_contract audit:partition_cmake PROPERTIES LABELS "unit;audit" TIMEOUT 120)

add_executable(audit_neutral_abi neutral_abi.c neutral_abi_fixture.f90)
target_include_directories(audit_neutral_abi PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(audit_neutral_abi PRIVATE gemini3d_mpi_c gemini3d_c gemini3d_mpi gemini3d MUMPS::MUMPS)
set_property(TARGET audit_neutral_abi PROPERTY LINKER_LANGUAGE CXX)
foreach(version IN ITEMS 0 21)
  if(version EQUAL 21 AND NOT gemini3d_msis2)
    continue()
  endif()
  add_test(NAME audit:neutral_abi_msis${version} COMMAND audit_neutral_abi ${version})
  set_tests_properties(audit:neutral_abi_msis${version} PROPERTIES
    LABELS "unit;audit" WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
endforeach()

# Make every public library declaration a linker dependency, not only the main program's subset.
file(STRINGS ${PROJECT_SOURCE_DIR}/include/gemini3d.h declarations REGEX "^extern void ")
set(link_source "#include \"gemini3d.h\"\n")
foreach(declaration IN LISTS declarations)
  string(REGEX REPLACE "^extern void ([A-Za-z0-9_]+).*" "\\1" symbol "${declaration}")
  string(APPEND link_source "decltype(&${symbol}) volatile link_${symbol} = &${symbol};\n")
endforeach()
string(APPEND link_source "int main() { return 0; }\n")
file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/public_symbols.cpp CONTENT "${link_source}")
add_executable(audit_public_symbols ${CMAKE_CURRENT_BINARY_DIR}/public_symbols.cpp)
target_include_directories(audit_public_symbols PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(audit_public_symbols PRIVATE gemini3d_mpi_c gemini3d_c gemini3d_mpi gemini3d MUMPS::MUMPS)
add_test(NAME audit:public_symbols COMMAND audit_public_symbols)
set_tests_properties(audit:public_symbols PROPERTIES LABELS "unit;audit")

if(Python_Interpreter_FOUND)
  add_test(NAME audit:filename_cadence COMMAND ${Python_EXECUTABLE}
    ${CMAKE_CURRENT_LIST_DIR}/test_filename_cadence.py --exe $<TARGET_FILE:audit_config>
    --work ${CMAKE_CURRENT_BINARY_DIR}/cadence-contract)
  set_tests_properties(audit:filename_cadence PROPERTIES LABELS "unit;audit" TIMEOUT 60)
endif()
if(NUMPY_FOUND AND H5PY_FOUND)
  add_test(NAME audit:launcher_preserves_outputs COMMAND ${Python_EXECUTABLE}
    ${CMAKE_CURRENT_LIST_DIR}/test_launcher_contract.py --exe $<TARGET_FILE:gemini3d.run>
    --work ${CMAKE_CURRENT_BINARY_DIR}/launcher-contract)
  set_tests_properties(audit:launcher_preserves_outputs PROPERTIES LABELS "unit;audit" TIMEOUT 60)
endif()
