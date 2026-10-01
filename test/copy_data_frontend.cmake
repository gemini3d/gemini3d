if(NOT IS_DIRECTORY ${outdir}/inputs)
  file(MAKE_DIRECTORY ${outdir})
  file(COPY ${refdir}/inputs DESTINATION ${outdir})
endif()
