# download reference data like:
#
# cmake -Dnames="2dns_fang;2dew_fang;3d_fang" -P cmake/download.cmake

function(gemini_download_ref_data names)

if(CMAKE_VERSION VERSION_LESS 3.19)
  file(READ ${CMAKE_CURRENT_LIST_DIR}/gemini3d_url.txt _refj)
else()
  file(READ ${CMAKE_CURRENT_LIST_DIR}/gemini3d_url.json _refj)
endif()

foreach(name ${names})

  if(CMAKE_VERSION VERSION_LESS 3.19)
    string(REGEX MATCH "${name} ([^ \t\r\n]*)" m ${_refj})
    set(url ${CMAKE_MATCH_1})
  else()
    string(JSON url GET ${_refj} ${name} url)
    # optional checksum
    string(JSON md5 ERROR_VARIABLE e GET ${_refj} ${name} md5)
  endif()

  set(archive_name test${name}.zip)
  set(ref_root ${CMAKE_CURRENT_LIST_DIR}/../tests/data)
  set(ref_dir ${ref_root}/test${name})
  set(archive ${ref_root}/${archive_name})

  if(EXISTS ${ref_dir}/inputs/config.nml)
    continue()
  endif()

  if(NOT EXISTS ${archive})
    if(md5)
      message(STATUS "download ${archive}  MD5: ${md5}")
      file(DOWNLOAD ${url} ${archive} TLS_VERIFY ON EXPECTED_HASH MD5=${md5})
    else()
      message(STATUS "download ${archive}")
      file(DOWNLOAD ${url} ${archive} TLS_VERIFY ON)
    endif()
  endif()

  message(STATUS "extract ref data to ${ref_dir}")
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.18)
    file(ARCHIVE_EXTRACT INPUT ${archive} DESTINATION ${ref_root})
  else()
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf ${archive} WORKING_DIRECTORY ${ref_root})
  endif()

endforeach()

endfunction(gemini_download_ref_data)

if(CMAKE_SCRIPT_MODE_FILE AND NOT CMAKE_PARENT_LIST_FILE)
  if(NOT names)
    message(FATAL_ERROR "specify reference data to download like:
      cmake -Dnames=\"2dns_fang;2dew_fang;3d_fang\" -P cmake/download.cmake")
  endif()

  gemini_download_ref_data("${names}")  # must have quotes to pass more than one value
endif()
