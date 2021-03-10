cmake_minimum_required(VERSION 3.13...3.20)

function(download_archive url archive)

if(md5)
  message(STATUS "download ${archive}  MD5: ${md5}")
  file(DOWNLOAD ${url} ${archive} TLS_VERIFY ON EXPECTED_HASH MD5=${md5})
else()
  message(STATUS "download ${archive}")
  file(DOWNLOAD ${url} ${archive} TLS_VERIFY ON)
endif()

endfunction(download_archive)


function(gemini_download_ref_data name)

# sanity check to avoid making mess
if(NOT refroot)
  message(FATAL_ERROR "must provide 'refroot' e.g. \${PROJECT_SOURCE_DIR}/test_data")
endif()

if(CMAKE_VERSION VERSION_LESS 3.19)
  file(READ ${CMAKE_CURRENT_LIST_DIR}/gemini3d_url.txt _refj)
else()
  file(READ ${CMAKE_CURRENT_LIST_DIR}/gemini3d_url.json _refj)
endif()

if(CMAKE_VERSION VERSION_LESS 3.19)
  string(REGEX MATCH "${name} ([^ \t\r\n]*)" m ${_refj})
  set(url ${CMAKE_MATCH_1})
else()
  string(JSON url GET ${_refj} ${name} url)
  # optional checksum
  string(JSON md5 ERROR_VARIABLE e GET ${_refj} ${name} md5)
endif()

set(archive_name test${name}.zip)
set(ref_dir ${refroot}/test${name})
set(archive ${refroot}/${archive_name})

# check if extracted data exists and is up to date
if(CMAKE_VERSION VERSION_LESS 3.19)
  # trivial check
  if(IS_DIRECTORY ${ref_dir})
    return()
  endif()
elseif(EXISTS ${ref_dir}/md5sum.txt AND md5)
  file(STRINGS ${ref_dir}/md5sum.txt _md5 REGEX "[a-f0-9]" LIMIT_INPUT 32 LENGTH_MAXIMUM 32 LIMIT_COUNT 1)

  if(_md5 STREQUAL ${md5})
    return()
  endif()
elseif(IS_DIRECTORY ${ref_dir})
  # missing md5, do trivial check
  return()
endif()

# check if archive .zip up to date
if(NOT EXISTS ${archive})
  download_archive(${url} ${archive})
endif()
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.19)
  file(MD5 ${archive} _md5)
  if(NOT _md5 STREQUAL ${md5})
    download_archive(${url} ${archive})
  endif()
endif()

message(STATUS "extract ref data to ${ref_dir}")
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.18)
  file(ARCHIVE_EXTRACT INPUT ${archive} DESTINATION ${refroot})
else()
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf ${archive} WORKING_DIRECTORY ${refroot})
endif()

file(MD5 ${archive} _md5)
file(WRITE ${ref_dir}/md5sum.txt ${_md5})

endfunction(gemini_download_ref_data)


gemini_download_ref_data(${name})

if(DEFINED outdir)
  # copy sim inputs into build/test${name}/inputs
  file(COPY ${refroot}/test${name}/inputs DESTINATION ${outdir})
endif()
