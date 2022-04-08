# this file defines simulation tests.
# The names of these must match those in ref_data.json

set(sim_timeout 1800)  # seconds to allow simulation tests to run

# --- setup tests
# NOTE: cpp just 2d is only to save test time. Works fine with 3d also.

set(gemini_test_sim_names)

if(hwm14)
  list(APPEND gemini_test_sim_names mini2dns_hwm14_fang mini2dew_hwm14_fang mini3d_hwm14_fang)
  if(cpp)
    list(APPEND gemini_test_sim_names mini2dns_hwm14_fang_cpp)
  endif()

  if(glow)
    list(APPEND gemini_test_sim_names mini2dns_hwm14_glow mini2dew_hwm14_glow mini3d_hwm14_glow)
    if(cpp)
      list(APPEND gemini_test_sim_names mini2dns_hwm14_glow_cpp)
    endif()
  endif(glow)
else()
  list(APPEND gemini_test_sim_names mini2dns_fang mini2dew_fang mini3d_fang)
  if(cpp)
    list(APPEND gemini_test_sim_names mini2dns_fang_cpp)
  endif()

  if(glow)
    list(APPEND gemini_test_sim_names mini2dns_glow mini2dew_glow mini3d_glow)
    if(cpp)
      list(APPEND gemini_test_sim_names mini2dns_glow_cpp)
    endif()
  endif(glow)
endif()


foreach(_s ${gemini_test_sim_names})
  setup_gemini_test(${_s} ${sim_timeout})

  if(python)
    setup_magcalc_test(${_s})
  endif()
endforeach()
