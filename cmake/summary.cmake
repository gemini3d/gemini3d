include(FeatureSummary)

# --- recommendations

if(CMAKE_GENERATOR MATCHES "Visual Studio")
  message(WARNING "Visual Studio generator ${CMAKE_GENERATOR} is not supported. Please use \"MinGW Makefiles\" or Ninja:

  cmake -G Ninja -B ${PROJECT_BINARY_DIR}
  ")
endif()

if(CMAKE_GENERATOR MATCHES "Ninja" AND CMAKE_VERSION VERSION_GREATER_EQUAL 3.27.0 AND CMAKE_VERSION VERSION_LESS 3.27.9)
  message(WARNING "CMake 3.27.0..3.27.8 has a bug with Ninja causing build failures.
  Suggest using CMake outside this range or:
  cmake -Bbuild -G \"Unix Makefiles\"
  ")
endif()

# --- options

add_feature_info(DevMode dev "Gemini developer mode")

add_feature_info(GLOW glow "airglow / aurora model")
add_feature_info(HWM14 hwm14 "HWM14 neutral winds model")

add_feature_info(PyGemini python "simulation generation, HPC script generator and plotting")
add_feature_info(MatGemini matlab "checks not as extensive as Python, and slow")

# print to screen
feature_summary(WHAT ENABLED_FEATURES DISABLED_FEATURES)
