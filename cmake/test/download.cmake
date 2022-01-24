cmake_minimum_required(VERSION 3.20...3.22)

include(${CMAKE_CURRENT_LIST_DIR}/../CheckTLS.cmake)
check_tls()

function(download_archive url archive exp_hash)

message(STATUS "DOWNLOAD: ${url} => ${archive}  sha256: ${exp_hash}")
file(DOWNLOAD ${url} ${archive} INACTIVITY_TIMEOUT 15)
file(SHA256 ${archive} hash)

if(hash STREQUAL ${exp_hash})
  return()
endif()

if(hash STREQUAL "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855")
  message(FATAL_ERROR "${url} failed to download: ${archive} is an empty file.")
endif()

message(FATAL_ERROR "${url} failed to download ${archive}")

endfunction(download_archive)


function(gemini_download_ref_data name refroot)

file(READ ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/ref_data.json _refj)

string(JSON url GET ${_refj} tests ${name} url)
string(JSON archive_name GET ${_refj} tests ${name} archive)
string(JSON hash GET ${_refj} tests ${name} sha256)

cmake_path(APPEND archive ${refroot} ${archive_name})

# check if extracted data exists and is up to date
if(EXISTS ${ref_dir}/sha256sum.txt)
  file(STRINGS ${ref_dir}/sha256sum.txt _hash REGEX "[a-f0-9]" LIMIT_INPUT 64 LENGTH_MAXIMUM 64 LIMIT_COUNT 1)

  if(_hash STREQUAL ${hash})
    return()
  else()
    message(STATUS "${name}: hash mismatch: expected ${hash} != ${_hash}")
  endif()
else()
  message(STATUS "${name}: missing hash file, seeing if we need to download and/or extract")
endif()

# check if archive up to date
if(NOT EXISTS ${archive})
  download_archive(${url} ${archive} ${hash})
endif()
file(SHA256 ${archive} _hash)
if(NOT _hash STREQUAL ${hash})
  download_archive(${url} ${archive} ${hash})
endif()

message(STATUS "EXTRACT: ${name}: ${archive} => ${ref_dir}")
file(ARCHIVE_EXTRACT INPUT ${archive} DESTINATION ${ref_dir})

file(SHA256 ${archive} _hash)
file(WRITE ${ref_dir}/sha256sum.txt ${_hash})

endfunction(gemini_download_ref_data)

# scripted part, needs to be in this order

cmake_path(APPEND ref_dir ${refroot} ${name})

gemini_download_ref_data(${name} ${refroot})

# copy sim inputs into build/${name}/inputs
file(COPY ${ref_dir}/inputs DESTINATION ${outdir})
