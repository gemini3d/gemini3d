# Test parameters (change for each test)
set(TESTDIR test2d)
set(REFNAME zenodo2d)
set(REFDIR ../simulations)
set(zenodoHash 5bd1bce1a465ccec5af813f8b7959ec8)
set(zenodoNumber 2394322)
set(firstfile 20130220_18000.000001.dat)
# --- ensure reference data is available for self-test
download_testfiles(${zenodoHash} ${zenodoNumber} ${REFNAME} ${PROJECT_SOURCE_DIR}/${REFDIR})

# --- test main exe
run_gemini_test(Gemini2D ${TESTDIR} 900)

# --- evaluate output accuracy vs. reference from Matt's HPC
compare_gemini_output(Compare2D ${TESTDIR} ${REFDIR}/${REFNAME} ${firstfile})

