# this file defines simulation tests.
# The names of these must match those in cmake/config/gemini3d_url.json

set(sim_timeout 1800)  # seconds to allow simulation tests to run

# --- setup tests

set(_tests 2dns_fang 2dew_fang 3d_fang)

if(glow)
  list(APPEND _tests 2dns_glow 2dew_glow 3d_glow)
endif(glow)

foreach(_s ${_tests})
  setup_gemini_test(${_s} ${sim_timeout})
endforeach()
