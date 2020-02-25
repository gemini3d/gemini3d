#== Fang
setup_gemini_test(
gemini:2d_fang
gemini.bin
test2d_fang
${CMAKE_CURRENT_SOURCE_DIR}/tests/data/test2d_fang
300)

compare_gemini_output(
gemini:compare:2d_fang
test2d_fang
${CMAKE_CURRENT_SOURCE_DIR}/tests/data/test2d_fang)

#== glow
if(glow)

setup_gemini_test(
gemini:2d_glow
gemini.bin
test2d_glow
${CMAKE_CURRENT_SOURCE_DIR}/tests/data/test2d_glow
300)

compare_gemini_output(
gemini:compare:2d_glow
test2d_glow
${CMAKE_CURRENT_SOURCE_DIR}/tests/data/test2d_glow)

endif(glow)
