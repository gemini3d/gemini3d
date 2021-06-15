# this file defines simulation tests.
# The names of these must match those in cmake/test/ref_data.json

set(sim_timeout 1800)  # seconds to allow simulation tests to run

# --- setup tests

set(_tests mini2dns_fang mini2dew_fang mini3d_fang)

if(glow)
  list(APPEND _tests mini2dns_glow mini2dew_glow mini3d_glow)
endif(glow)

foreach(_s ${_tests})
  setup_gemini_test(${_s} ${sim_timeout})

  if(PYGEMINI_DIR)
    setup_magcalc_test(${_s})
  endif()
endforeach()
