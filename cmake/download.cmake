function(download_testfiles HASH REFNUM REFNAME ROOT)
# FetchContent is too aggressive, it deletes the output directory before extracting.
# That could delete wrong directory causing data loss.
# So use this more specific but safe method.

set(ARCHIVE ${REFNAME}.zip)
set(URL https://zenodo.org/record/${REFNUM}/files/${REFNAME}.zip?download=1)

# --- ensure reference data is available for self-test
if(EXISTS ${ROOT}/${REFNAME})
  return()
endif()

if(NOT EXISTS ${ROOT}/${ARCHIVE})

file(DOWNLOAD ${URL} ${ROOT}/${ARCHIVE})
# don't check hash during download in case of failure, so it doesn't stop using rest of program.
endif()

file(MD5 ${ROOT}/${ARCHIVE} FHASH)

if(FHASH STREQUAL HASH)
# extract zip file (os-agnostic)
execute_process(COMMAND ${CMAKE_COMMAND} -E tar -xf ${ARCHIVE}
                WORKING_DIRECTORY ${ROOT})
else()
  message(WARNING "Self-tests may not work. Reference data did not match hash ${HASH} \n Download ${URL}  to  ${ARCHIVE} \n")
endif()


endfunction(download_testfiles)
