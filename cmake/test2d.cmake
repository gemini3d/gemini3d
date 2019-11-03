download_testfiles(
57d72fd0005247c8eedf122ac4670ad0
3477385
zenodo2d_fang
${PROJECT_SOURCE_DIR}/tests/data
)

setup_gemini_test(
Gemini2d_fang
gemini_fang.bin
test2d_fang
tests/data/zenodo2d_fang
300
)

compare_gemini_output(
Compare2d_fang
test2d_fang
tests/data/zenodo2d_fang
20130220_18000.000001
)

if(useglow)
download_testfiles(
557bc6a91d8bf3464abdc5c8784f3042
3477385
zenodo2d_glow
${PROJECT_SOURCE_DIR}/tests/data
)

setup_gemini_test(
Gemini2d_glow
gemini_glow.bin
test2d_glow
tests/data/zenodo2d_glow
300
)

compare_gemini_output(
Compare2d_glow
test2d_glow
tests/data/zenodo2d_glow
20130220_18000.000001
)
endif(useglow)
