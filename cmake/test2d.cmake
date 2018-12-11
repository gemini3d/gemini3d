# Test parameters (change for each test)
set(TESTDIR test2d)
set(REFNAME zenodo2d)
set(REFDIR ../simulations)
set(zenodoHash c577e5b9a9b6e0d8506b1798181cf230)
set(zenodoNumber 1975078)
# --- ensure reference data is available for self-test
download_testfiles(${zenodoHash}
                   ${zenodoNumber}
                   ${REFNAME}
                   ${PROJECT_SOURCE_DIR}/${REFDIR})

# --- test main exe
run_gemini_test(Gemini2D ${TESTDIR} 900)

# --- evaluate output accuracy vs. reference from Matt's HPC
compare_gemini_output(Compare2D ${TESTDIR} ${REFDIR}/${REFNAME})

