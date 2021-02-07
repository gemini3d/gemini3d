include(${CMAKE_CURRENT_LIST_DIR}/download.cmake)

if(IS_DIRECTORY ${outdir})
  return()
endif()

if(NOT IS_DIRECTORY ${refdir})
  gemini_download_ref_data("${testname}")
endif()

file(COPY ${refdir} DESTINATION ${outdir})
