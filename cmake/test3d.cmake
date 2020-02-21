#== Fang
download_testfiles(
todo
https://zenodo.invalid/record/
test3d_fang
${CMAKE_SOURCE_DIR}/tests/data
)

setup_gemini_test(
Gemini3d_fang
gemini.bin
test3d_fang
${CMAKE_SOURCE_DIR}/tests/data/test3d_fang
600
)

compare_gemini_output(
Compare3d_fang
test3d_fang
${CMAKE_SOURCE_DIR}/tests/data/test3d_fang
)

#== Glow
if(glow)
download_testfiles(
todo
https://zenodo.invalid/record/
test3d_glow
${CMAKE_SOURCE_DIR}/tests/data
)

setup_gemini_test(
Gemini3d_glow
gemini.bin
test3d_glow
${CMAKE_SOURCE_DIR}/tests/data/test3d_glow
1800
)

compare_gemini_output(
Compare3d_glow
test3d_glow
${CMAKE_SOURCE_DIR}/tests/data/test3d_glow
)
endif(glow)
