#== Fang
download_testfiles(
57d72fd0005247c8eedf122ac4670ad0
3477385
zenodo2d_fang
${CMAKE_SOURCE_DIR}/tests/data
)

setup_gemini_test(
Gemini2d_fang
gemini_fang.bin
test2d_fang
${CMAKE_SOURCE_DIR}/tests/data/zenodo2d_fang
300
)

compare_gemini_output(
Compare2d_fang
test2d_fang
${CMAKE_SOURCE_DIR}/tests/data/zenodo2d_fang
)

#== glow
if(useglow)
download_testfiles(
557bc6a91d8bf3464abdc5c8784f3042
3477385
zenodo2d_glow
${CMAKE_SOURCE_DIR}/tests/data
)

setup_gemini_test(
Gemini2d_glow
gemini_glow.bin
test2d_glow
${CMAKE_SOURCE_DIR}/tests/data/zenodo2d_glow
300
)

compare_gemini_output(
Compare2d_glow
test2d_glow
${CMAKE_SOURCE_DIR}/tests/data/zenodo2d_glow
)
endif(useglow)
