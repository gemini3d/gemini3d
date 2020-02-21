function(download_testfiles HASH URL REFNAME ROOT)
# FetchContent is too aggressive, it deletes the output directory before extracting.
# That could delete wrong directory causing data loss.
# So use this more specific but safe method.

# We can't use CMAKE_ROLE to defer downloading because the tests parameters are defined
# using the test reference data.

set(_zipname ${REFNAME}.zip)
set(_zipfile ${ROOT}/${_zipname})

if(EXISTS ${ROOT}/${REFNAME})
  # reference file already extracted to directory
  return()
endif()

if(NOT EXISTS ${_zipfile})
  # don't check hash during download in case of failure, so it doesn't stop using rest of program.
  file(DOWNLOAD ${URL} ${_zipfile})
endif()

file(SIZE ${_zipfile} _zipsize)
if(_zipsize LESS 1000000)
  message(STATUS "SKIP: outdated URL or truncated download: ${REFNAME}")
  return()
endif()

file(MD5 ${_zipfile} _ziphash)

if(NOT _ziphash STREQUAL HASH) # os-agnostic
  message(STATUS "SKIP: wrong URL or corrupted file: ${REFNAME}")
  return()
endif()

execute_process(COMMAND ${CMAKE_COMMAND} -E tar -xf ${_zipname} WORKING_DIRECTORY ${ROOT})


endfunction(download_testfiles)
