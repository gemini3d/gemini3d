include(FeatureSummary)

if(gemini3d_IS_TOP_LEVEL)

add_feature_info(GLOW gemini3d_glow "airglow / aurora model")
add_feature_info(HWM14 gemini3d_hwm14 "HWM14 neutral winds model")

add_feature_info(PyGemini gemini3d_python "simulation generation, HPC script generator and plotting")
add_feature_info(MatGemini gemini3d_matlab "checks not as extensive as Python, and slow")

# print to screen
feature_summary(WHAT ENABLED_FEATURES DISABLED_FEATURES)

endif()
