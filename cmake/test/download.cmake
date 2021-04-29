cmake_minimum_required(VERSION 3.13...3.20)

function(download_archive url archive)

if(sha256)
  message(STATUS "download ${archive}  sha256: ${sha256}")
  file(DOWNLOAD ${url} ${archive} EXPECTED_HASH SHA256=${sha256})
else()
  message(STATUS "download ${archive}")
  file(DOWNLOAD ${url} ${archive})
endif()

endfunction(download_archive)


function(gemini_download_ref_data name)

# sanity check to avoid making mess
if(NOT refroot)
  message(FATAL_ERROR "must provide 'refroot' e.g. \${PROJECT_SOURCE_DIR}/test_data")
endif()

file(READ ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/ref_data.json _refj)

string(JSON url GET ${_refj} tests ${name} url)
string(JSON archive_name GET ${_refj} tests ${name} archive)
# optional checksum
string(JSON sha256 ERROR_VARIABLE e GET ${_refj} tests ${name} sha256)

cmake_path(APPEND archive ${refroot} ${archive_name})

# check if extracted data exists and is up to date
if(EXISTS ${ref_dir}/sha256sum.txt AND sha256)
  file(STRINGS ${ref_dir}/sha256sum.txt _sha256 REGEX "[a-f0-9]" LIMIT_INPUT 32 LENGTH_MAXIMUM 32 LIMIT_COUNT 1)

  if(_sha256 STREQUAL ${sha256})
    return()
  endif()
elseif(IS_DIRECTORY ${ref_dir})
  # missing has, do trivial check
  return()
endif()

# check if archive up to date
if(NOT EXISTS ${archive})
  download_archive(${url} ${archive})
endif()
file(SHA256 ${archive} _sha256)
if(NOT _sha256 STREQUAL ${sha256})
  download_archive(${url} ${archive})
endif()

message(STATUS "extract ref data to ${ref_dir}")
file(ARCHIVE_EXTRACT INPUT ${archive} DESTINATION ${ref_dir})

file(SHA256 ${archive} _sha256)
file(WRITE ${ref_dir}/sha256sum.txt ${_sha256})

endfunction(gemini_download_ref_data)

# scripted part, needs to be in this order

cmake_path(APPEND ref_dir ${refroot} ${name})

gemini_download_ref_data(${name})

# copy sim inputs into build/${name}/inputs
message(STATUS "TRACE: ${ref_dir} ${outdir}")
file(COPY ${ref_dir}/inputs DESTINATION ${outdir})
