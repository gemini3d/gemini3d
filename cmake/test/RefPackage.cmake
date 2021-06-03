function(ref_package data_dir ref_json_file name)

set(ARC_TYPE zst)

cmake_path(GET ref_json_file PARENT_PATH arc_dir)
cmake_path(APPEND archive ${arc_dir} ${name}.${ARC_TYPE})

add_test(NAME "package:archive:${name}"
  COMMAND ${CMAKE_COMMAND} -Din:PATH=${data_dir} -Dout:FILEPATH=${archive} -Dref_json_file:FILEPATH=${ref_json_file} -Dname=${name} -P ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/archive.cmake)

set_tests_properties("package:archive:${name}" PROPERTIES
FIXTURES_REQUIRED hdf5:${name}
# FIXTURES_SETUP ${name}:upload_fxt
LABELS package
REQUIRED_FILES "${data_dir}/inputs/config.nml;${data_dir}/output.nml"
RESOURCE_LOCK archive_lock  # else will race overwrite ref.json
TIMEOUT 30)

endfunction(ref_package)
