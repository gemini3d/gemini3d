# Test parameters (change for each test)
set(TESTDIR test2d_glow)
set(REFNAME zenodo2d_glow)
set(REFDIR ../simulations)
set(zenodoHash c5bbbbff3bdde85b6d7e9470bc3751a2)
set(zenodoNumber 2520780)
set(firstfile 20130220_18000.000001.dat)
# --- ensure reference data is available for self-test
download_testfiles(${zenodoHash} ${zenodoNumber} ${REFNAME} ${PROJECT_SOURCE_DIR}/${REFDIR})

setup_gemini_test(GeminiGlow2D ${TESTDIR} ${REFDIR}/${REFNAME} 300)

compare_gemini_output(CompareGlow2D ${TESTDIR} ${REFDIR}/${REFNAME} ${firstfile})
