download_testfiles(
225853d43937a70c9ef6726f90666645
2520920
zenodo3d
${PROJECT_SOURCE_DIR}/tests/data
)

setup_gemini_test(
Gemini3d_fang
gemini_fang.bin
test3d
tests/data/zenodo3d
600
)

compare_gemini_output(
Compare3d_fang
test3d
tests/data/zenodo3d
20130220_18000.000001.dat
)

if(useglow)
download_testfiles(
2c54bfde8aff0fb72d61115f04c361a7
2520920
zenodo3d_glow
${PROJECT_SOURCE_DIR}/tests/data
)

setup_gemini_test(
Gemini3d_glow
gemini_glow.bin
test3d_glow
tests/data/zenodo3d_glow
1200
)

compare_gemini_output(
Compare3d_glow
test3d_glow
tests/data/zenodo3d_glow
20130220_18000.000001.dat
)
endif(useglow)
