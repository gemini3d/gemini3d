# Test parameters (change for each test)
set(TESTDIR test3d)
set(REFNAME zenodo3d)
set(REFDIR ../simulations)
set(zenodoHash 225853d43937a70c9ef6726f90666645)
set(zenodoNumber 1473195)
set(firstfile 20130220_18000.000001.dat)
# --- ensure reference data is available for self-test
download_testfiles(${zenodoHash} ${zenodoNumber} ${REFNAME} ${PROJECT_SOURCE_DIR}/${REFDIR})

# --- test GEMINI
run_gemini_test(Gemini3D ${TESTDIR} 900)

# --- test GEMINI + GLOW
run_gemini_test(GlowGemini3D ${TESTDIR}_glow 900)
    
# --- evaluate output accuracy vs. reference from Matt's HPC
compare_gemini_output(Compare3D ${TESTDIR} ${REFDIR}/${REFNAME} ${firstfile})

