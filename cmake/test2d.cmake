#== Fang
download_testfiles(
03c183bbc91706223313e5c15771918e
https://zenodo.org/record/3677638/files/test2d_fang.zip?download=1
test2d_fang
${CMAKE_SOURCE_DIR}/tests/data
)

setup_gemini_test(
gemini:2d_fang
gemini.bin
test2d_fang
${CMAKE_SOURCE_DIR}/tests/data/test2d_fang
300
)

compare_gemini_output(
gemini:compare:2d_fang
test2d_fang
${CMAKE_SOURCE_DIR}/tests/data/test2d_fang
)

#== glow
if(glow)
download_testfiles(
bd9a9c38bb462cc22cc0ea232e03dc21
https://zenodo.org/record/3677638/files/test2d_glow.zip?download=1
test2d_glow
${CMAKE_SOURCE_DIR}/tests/data
)

setup_gemini_test(
gemini:2d_glow
gemini.bin
test2d_glow
${CMAKE_SOURCE_DIR}/tests/data/test2d_glow
300
)

compare_gemini_output(
gemini:compare:2d_glow
test2d_glow
${CMAKE_SOURCE_DIR}/tests/data/test2d_glow
)
endif(glow)
