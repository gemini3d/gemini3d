# this file defines simulation tests.
# The names of these must match those in pygemini/src/gemini3d/tests/gemini3d_url.ini
# we are trying to avoid multiple copies of that URL/MD5 data.

# === Fang

# North-South 2D
setup_gemini_test(2dns_fang 300)

# East-West 2D
setup_gemini_test(2dew_fang 300)

# 3D
setup_gemini_test(3d_fang 600)

# === glow

if(glow)

  # North-South 2D
  setup_gemini_test(2dns_glow 300)

  # East-West 2D
  setup_gemini_test(2dew_glow 300)

  # 3D
  setup_gemini_test(3d_glow 1800)

endif(glow)
