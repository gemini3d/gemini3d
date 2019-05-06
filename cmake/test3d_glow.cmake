# Test parameters (change for each test)
set(TESTDIR test3d_glow)
set(REFNAME zenodo3d_glow)
set(REFDIR ../simulations)
set(zenodoHash 2c54bfde8aff0fb72d61115f04c361a7)
set(zenodoNumber 2520920)
set(firstfile 20130220_18000.000001.dat)
# --- ensure reference data is available for self-test
download_testfiles(${zenodoHash} ${zenodoNumber} ${REFNAME} ${PROJECT_SOURCE_DIR}/${REFDIR})

setup_gemini_test(GeminiGlow3D ${TESTDIR} ${REFDIR}/${REFNAME} 1200)

compare_gemini_output(CompareGlow3D ${TESTDIR} ${REFDIR}/${REFNAME} ${firstfile})
