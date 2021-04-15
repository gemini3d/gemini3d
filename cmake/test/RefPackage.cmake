function(ref_package data_dir ref_json_file name)

set(ARC_TYPE zstd)

cmake_path(GET ref_json_file PARENT_PATH arc_dir)
cmake_path(APPEND archive ${arc_dir} ${name}.${ARC_TYPE})

if(NOT EXISTS ${ref_json_file})
  # we will read this file and overwrite with new JSON for each sim
  file(MAKE_DIRECTORY ${arc_dir})
  file(WRITE ${ref_json_file} "{}")
endif()

add_test(NAME "package:archive:${name}"
  COMMAND ${CMAKE_COMMAND} -Din:PATH=${data_dir} -Dout:FILEPATH=${archive} -Dref_json_file:FILEPATH=${ref_json_file} -Dname=${name} -P ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/archive.cmake)

set_tests_properties("package:archive:${name}" PROPERTIES
FIXTURES_REQUIRED hdf5:${name}
# FIXTURES_SETUP ${name}:upload_fxt
LABELS package
REQUIRED_FILES "${data_dir}/inputs/config.nml;${data_dir}/output.nml"
TIMEOUT 30)

endfunction(ref_package)
