include(FeatureSummary)

if(gemini3d_IS_TOP_LEVEL)
add_feature_info(DevMode dev "Gemini developer mode")

add_feature_info(GLOW glow "airglow / aurora model")
add_feature_info(HWM14 hwm14 "HWM14 neutral winds model")

add_feature_info(PyGemini python "simulation generation, HPC script generator and plotting")
add_feature_info(MatGemini matlab "checks not as extensive as Python, and slow")

# print to screen
feature_summary(WHAT ENABLED_FEATURES DISABLED_FEATURES)
endif()
