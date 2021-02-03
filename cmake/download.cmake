# download reference data like:
#
# cmake -Ddownload="2dns_fang;2dew_fang;3d_fang" -P cmake/download.cmake

if(CMAKE_VERSION VERSION_LESS 3.19)
  file(READ ${CMAKE_CURRENT_LIST_DIR}/gemini3d_url.txt _refj)
else()
  file(READ ${CMAKE_CURRENT_LIST_DIR}/gemini3d_url.json _refj)
endif()


foreach(d ${download})

  if(CMAKE_VERSION VERSION_LESS 3.19)
    string(REGEX MATCH "${d} ([^ \t\r\n]*)" m ${_refj})
    set(url ${CMAKE_MATCH_1})
  else()
    string(JSON url GET ${_refj} ${d} url)
  endif()

  set(archive_name test${d}.zip)
  set(ref_root ${CMAKE_CURRENT_LIST_DIR}/../tests/data)
  set(ref_dir ${ref_root}/test${d})
  set(archive ${ref_root}/${archive_name})

  if(EXISTS ${ref_dir}/inputs/config.nml)
    continue()
  endif()

  if(NOT EXISTS ${archive})
    message(STATUS "download ${archive}")

    file(DOWNLOAD ${url} ${archive} TLS_VERIFY ON)
  endif()

  message(STATUS "extract ref data to ${ref_dir}")
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.18)
    file(ARCHIVE_EXTRACT INPUT ${archive} DESTINATION ${ref_root})
  else()
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf ${archive} WORKING_DIRECTORY ${ref_root})
  endif()

endforeach()
