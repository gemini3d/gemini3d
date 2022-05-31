set_directory_properties(PROPERTIES LABELS io)


cmake_path(SET _test_file ${PROJECT_SOURCE_DIR}/test/config/config_example.nml)

add_executable(test_namelist test_nml.f90)
add_test(NAME Fortran_namelist COMMAND test_namelist ${_test_file})
set_tests_properties(Fortran_namelist PROPERTIES
REQUIRED_FILES ${_test_file}
LABELS unit
TIMEOUT 10
)

add_executable(test_compiler_vendor test_compiler_vendor.f90)
target_link_libraries(test_compiler_vendor PRIVATE config
ffilesystem::filesystem
)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  set(compiler_vendor Intel)
else()
  set(compiler_vendor ${CMAKE_Fortran_COMPILER_ID})
endif()

add_test(NAME compiler_vendor COMMAND test_compiler_vendor)
set_tests_properties(compiler_vendor PROPERTIES
PASS_REGULAR_EXPRESSION ${compiler_vendor}
LABELS unit
TIMEOUT 10
)

add_executable(test_expand_envvar test_expand_envvar.f90)
target_link_libraries(test_expand_envvar PRIVATE config
ffilesystem::filesystem
)

add_test(NAME expand_envvar COMMAND test_expand_envvar)
set_tests_properties(expand_envvar PROPERTIES
ENVIRONMENT test__gem=hello
LABELS unit
TIMEOUT 10
)

# --- test props

if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  set_target_properties(test_compiler_vendor test_expand_envvar PROPERTIES LINKER_LANGUAGE Fortran)
else()
  set_target_properties(test_compiler_vendor test_expand_envvar PROPERTIES LINKER_LANGUAGE CXX)
endif()

dll_test_path(ffilesystem::filesystem "expand_envvar;compiler_vendor")


if(f2018assumed_rank)
  add_library(debug_dump OBJECT debug_dump.f90)
  target_link_libraries(debug_dump PRIVATE const h5fortran::h5fortran)
  set_target_properties(debug_dump PROPERTIES EXCLUDE_FROM_ALL true)
endif()

add_executable(test_dump_write test_dump_write.f90)
target_link_libraries(test_dump_write PRIVATE errors const
h5fortran::h5fortran
$<$<BOOL:${HDF5_HAVE_PARALLEL}>:MPI::MPI_Fortran>
)

add_executable(test_dump_read test_dump_read.f90)
target_link_libraries(test_dump_read PRIVATE
h5fortran::h5fortran
$<$<BOOL:${HDF5_HAVE_PARALLEL}>:MPI::MPI_Fortran>
)

foreach(t root worker input)
  set(f ${CMAKE_CURRENT_BINARY_DIR}/dump_${t}.h5)

  add_test(NAME gemini:dump:${t} COMMAND test_dump_write ${t} ${f})
  set_tests_properties(gemini:dump:${t} PROPERTIES
  FIXTURES_SETUP gemini:dump_fxt
  WILL_FAIL true
  LABELS unit
  TIMEOUT 10
  )

  add_test(NAME gemini:dump:check:${t} COMMAND test_dump_read testing:${t} ${f})
  set_tests_properties(gemini:dump:check:${t} PROPERTIES
  FIXTURES_REQUIRED gemini:dump_fxt
  REQUIRED_FILES ${f}
  LABELS unit
  TIMEOUT 10
  )

  dll_test_path("h5fortran::h5fortran;HDF5::HDF5" "gemini:dump:${t};gemini:dump:check:${t}")
endforeach()


if(netcdf AND ${PROJECT_NAME}_BUILD_TESTING)
  add_executable(test_simsize_nc4 test_simsize_nc4.f90)
  target_link_libraries(test_simsize_nc4 PRIVATE nc4fortran::nc4fortran)
  # add_test(NAME simsize_nc4 COMMAND test_simsize_nc4)
  # set_tests_properties(simsize_nc4 PROPERTIES LABELS unit)
endif()
