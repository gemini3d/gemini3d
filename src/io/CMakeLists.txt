set_directory_properties(PROPERTIES LABELS io)

add_library(logging OBJECT logging.f90)

add_library(sanity_check OBJECT sanity_check.f90 $<TARGET_OBJECTS:errors>)
target_link_libraries(sanity_check PRIVATE const)

# --- error handling -- dump data to disk
add_library(errors OBJECT error_stop.f90)
target_link_libraries(errors PRIVATE const h5fortran::h5fortran)

if(${PROJECT_NAME}_BUILD_TESTING)

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

    add_test(NAME gemini:dump:${t} COMMAND $<TARGET_FILE:test_dump_write> ${t} ${f})
    set_tests_properties(gemini:dump:${t} PROPERTIES
    FIXTURES_SETUP gemini:dump_fxt
    WILL_FAIL true
    LABELS unit
    TIMEOUT 10
    )

    add_test(NAME gemini:dump:check:${t} COMMAND $<TARGET_FILE:test_dump_read> testing:${t} ${f})
    set_tests_properties(gemini:dump:check:${t} PROPERTIES
    FIXTURES_REQUIRED gemini:dump_fxt
    REQUIRED_FILES ${f}
    LABELS unit
    TIMEOUT 10
    )
  endforeach()

endif()

# --- config
add_library(config OBJECT config.f90 config_nml.f90 config_ini.f90)
target_link_libraries(config PRIVATE const
ffilesystem::filesystem ${lib_filesystem}
)

set(_test_file ${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests/config/config_example.nml)

if(${PROJECT_NAME}_BUILD_TESTING)

add_executable(test_namelist test_nml.f90)
add_test(NAME unit:namelist COMMAND $<TARGET_FILE:test_namelist> ${_test_file})
set_tests_properties(unit:namelist PROPERTIES
REQUIRED_FILES ${_test_file}
LABELS unit
TIMEOUT 10
)

add_executable(test_compiler_vendor test_compiler_vendor.f90)
target_link_libraries(test_compiler_vendor PRIVATE config
ffilesystem::filesystem ${lib_filesystem}
)
if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  set_target_properties(test_compiler_vendor PROPERTIES LINKER_LANGUAGE Fortran)
else()
  set_target_properties(test_compiler_vendor PROPERTIES LINKER_LANGUAGE CXX)
endif()
add_test(NAME unit:compiler_vendor
COMMAND $<TARGET_FILE:test_compiler_vendor>
)
set_tests_properties(unit:compiler_vendor PROPERTIES
PASS_REGULAR_EXPRESSION ${CMAKE_Fortran_COMPILER_ID}
LABELS unit
TIMEOUT 10
)

add_executable(test_expand_envvar test_expand_envvar.f90)
target_link_libraries(test_expand_envvar PRIVATE config
ffilesystem::filesystem ${lib_filesystem}
)
if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  set_target_properties(test_expand_envvar PROPERTIES LINKER_LANGUAGE Fortran)
else()
  set_target_properties(test_expand_envvar PROPERTIES LINKER_LANGUAGE CXX)
endif()
add_test(NAME unit:expand_envvar COMMAND $<TARGET_FILE:test_expand_envvar>)
set_tests_properties(unit:expand_envvar PROPERTIES
ENVIRONMENT test__gem=hello
LABELS unit
TIMEOUT 10
)

endif()

# --- io
configure_file(output.in.f90 output.f90 @ONLY)
add_library(io OBJECT io.f90 ${CMAKE_CURRENT_BINARY_DIR}/output.f90
aurora.f90 aurora_raw.f90 aurora_hdf5.f90 aurora_ncdf.f90
cond.f90 cond_hdf5.f90
mag.f90 mag_hdf5.f90 mag_ncdf.f90 mag_raw.f90
plasma_input.f90 plasma_input_hdf5.f90 plasma_input_ncdf.f90 plasma_input_raw.f90
plasma_output.f90 plasma_output_hdf5.f90 plasma_output_ncdf.f90 plasma_output_raw.f90
milestone.f90
)

target_link_libraries(io PRIVATE calculus config const reader grid mpimod sanity_check timeutils
h5fortran::h5fortran nc4fortran::nc4fortran
ffilesystem::filesystem ${lib_filesystem}
)

# --- reader
add_library(reader OBJECT reader.f90 reader_hdf5.f90 reader_ncdf.f90 reader_raw.f90)
target_link_libraries(reader PRIVATE const
h5fortran::h5fortran nc4fortran::nc4fortran
ffilesystem::filesystem ${lib_filesystem}
)

if(netcdf AND ${PROJECT_NAME}_BUILD_TESTING)
  add_executable(test_simsize_nc4 test_simsize_nc4.f90)
  target_link_libraries(test_simsize_nc4 PRIVATE nc4fortran::nc4fortran)
  # add_test(NAME unit:simsize_nc4 COMMAND test_simsize_nc4)
endif()
