function(make_archive in out ref_json_file name)

get_filename_component(in ${in} ABSOLUTE)
get_filename_component(out ${out} ABSOLUTE)
cmake_path(GET out EXTENSION LAST_ONLY ARC_TYPE)
cmake_path(GET out FILENAME archive_name)

if(ARC_TYPE STREQUAL .zst OR ARC_TYPE STREQUAL .zstd)
  # not usable due to internal paths always relative to PROJECT_BINARY_DIR
  # https://gitlab.kitware.com/cmake/cmake/-/issues/21653
  # file(ARCHIVE_CREATE
  #   OUTPUT ${out}
  #   PATHS ${in}
  #   COMPRESSION Zstd
  #   COMPRESSION_LEVEL 3)

  # need working_directory ${in} to avoid computer-specific relative paths
  # use . not ${in} as last argument to avoid more relative path issues
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar c ${out} --zstd .
    WORKING_DIRECTORY ${in})

elseif(ARC_TYPE STREQUAL .zip)

  execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar c ${out} --format=zip .
    WORKING_DIRECTORY ${in})

else()
  message(FATAL_ERROR "unknown archive type ${ARC_TYPE}")
endif()

# ensure a file was created (weak validation)
if(NOT EXISTS ${out})
  message(FATAL_ERROR "Archive ${out} was not created.")
endif()

file(SIZE ${out} fsize)
if(fsize LESS 10000)
  message(FATAL_ERROR "Archive ${out} may be malformed.")
endif()

# put hash in JSON
file(SHA256 ${out} sha256_hash)
file(READ ${ref_json_file} ref_json)
# will manually get URL after upload
string(JSON m ERROR_VARIABLE e GET ${ref_json} tests)

if(NOT m)
  string(JSON ref_json SET ${ref_json} tests "{}")
endif()
string(JSON ref_json SET ${ref_json} tests ${name} "{}")
string(JSON ref_json SET ${ref_json} tests ${name} url \"\")
string(JSON ref_json SET ${ref_json} tests ${name} archive \"${archive_name}\")
string(JSON ref_json SET ${ref_json} tests ${name} sha256 \"${sha256_hash}\")
file(WRITE ${ref_json_file} ${ref_json})

endfunction(make_archive)


make_archive(${in} ${out} ${ref_json_file} ${name})
