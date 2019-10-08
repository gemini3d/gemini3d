download_testfiles(
cf73d6eb166369c522da7a371492a1ce
3477330
zenodo3d_fang
${PROJECT_SOURCE_DIR}/tests/data
)

setup_gemini_test(
Gemini3d_fang
gemini_fang.bin
test3d_fang
tests/data/zenodo3d
600
)

compare_gemini_output(
Compare3d_fang
test3d_fang
tests/data/zenodo3d
20130220_18000.000001.dat
)

if(useglow)
download_testfiles(
3528946525295cc8271aa41bc262d7f1
3477330
zenodo3d_glow
${PROJECT_SOURCE_DIR}/tests/data
)

setup_gemini_test(
Gemini3d_glow
gemini_glow.bin
test3d_glow
tests/data/zenodo3d_glow
1800
)

compare_gemini_output(
Compare3d_glow
test3d_glow
tests/data/zenodo3d_glow
20130220_18000.000001.dat
)
endif(useglow)
