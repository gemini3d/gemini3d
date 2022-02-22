file(READ ${CMAKE_CURRENT_LIST_DIR}/libraries.json _libj)

if(matlab)
  string(JSON matgemini_git GET ${_libj} matgemini git)
endif()

# --- download reference data JSON file (for previously generated data)
cmake_path(APPEND arc_json_file ${CMAKE_CURRENT_BINARY_DIR} ref_data.json)
if(NOT EXISTS ${arc_json_file})
  string(JSON url GET ${_libj} ref_data url)
  file(DOWNLOAD ${url} ${arc_json_file} INACTIVITY_TIMEOUT 15)
endif()
