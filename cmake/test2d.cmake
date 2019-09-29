download_testfiles(
57d72fd0005247c8eedf122ac4670ad0
3464571
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
20130220_18000.000001.dat
)

if(useglow)
download_testfiles(
c5bbbbff3bdde85b6d7e9470bc3751a2
2520780
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
20130220_18000.000001.dat
)
endif(useglow)
