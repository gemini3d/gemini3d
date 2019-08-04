download_testfiles(
5bd1bce1a465ccec5af813f8b7959ec8
2520780
zenodo2d
${PROJECT_SOURCE_DIR}/tests/data
)

setup_gemini_test(
Gemini2d_fang
gemini_fang.bin
test2d
tests/data/zenodo2d
300
)

compare_gemini_output(
Compare2d_fang
test2d
tests/data/zenodo2d
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
