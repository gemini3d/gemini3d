# data is one-time downloaded using cmake -P cmake/download.cmake

#== Fang
setup_gemini_test(
gemini:3d_fang
gemini.bin
test3d_fang
${CMAKE_CURRENT_SOURCE_DIR}/tests/data/test3d_fang
600)

compare_gemini_output(
gemini:compare:3d_fang
test3d_fang
${CMAKE_CURRENT_SOURCE_DIR}/tests/data/test3d_fang)

#== Glow
if(glow)

setup_gemini_test(
gemini:3d_glow
gemini.bin
test3d_glow
${CMAKE_CURRENT_SOURCE_DIR}/tests/data/test3d_glow
1800)

compare_gemini_output(
gemini:compare:3d_glow
test3d_glow
${CMAKE_CURRENT_SOURCE_DIR}/tests/data/test3d_glow)

endif(glow)
