#== Fang
setup_gemini_test(2dew_fang 300)

setup_gemini_test(3d_fang 600)

# === glow
if(glow)
  setup_gemini_test(2dew_glow 300)

  setup_gemini_test(3d_glow 1800)
endif(glow)