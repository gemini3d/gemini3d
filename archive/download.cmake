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
  message(VERBOSE " ${ROOT}/${REFNAME} exists, skipping download")
  return()
endif()

if(NOT EXISTS ${_zipfile})
  # don't check hash during download in case of failure, so it doesn't stop using rest of program.
  file(DOWNLOAD ${URL} ${_zipfile}
    SHOW_PROGRESS)
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


# these are meant to be called as cmake -P
if(download_2d_fang)
download_testfiles(
  03c183bbc91706223313e5c15771918e
  https://zenodo.org/record/3677638/files/test2d_fang.zip?download=1
  test2d_fang
  ${CMAKE_CURRENT_SOURCE_DIR}/tests/data)
endif()

if(download_2d_glow)
download_testfiles(
  bd9a9c38bb462cc22cc0ea232e03dc21
  https://zenodo.org/record/3677638/files/test2d_glow.zip?download=1
  test2d_glow
  ${CMAKE_CURRENT_SOURCE_DIR}/tests/data
  )
endif()

if(download_3d_fang)
download_testfiles(
  b4c5fc43243b33549b8324c9a56ee198
  https://zenodo.org/record/3687202/files/test3d_fang.zip?download=1
  test3d_fang
  ${CMAKE_CURRENT_SOURCE_DIR}/tests/data)
endif()

if(download_3d_glow)
download_testfiles(
  d70c8ee5a699ae028b0ffecb750fb5c6
  https://zenodo.org/record/3687202/files/test3d_glow.zip?download=1
  test3d_glow
  ${CMAKE_CURRENT_SOURCE_DIR}/tests/data)
endif()
