include(FeatureSummary)

# --- recommendations

if(CMAKE_GENERATOR MATCHES "Visual Studio")
  message(WARNING "Visual Studio generator ${CMAKE_GENERATOR} is not supported. Please use \"MinGW Makefiles\" or Ninja:

  cmake -G Ninja -B ${PROJECT_BINARY_DIR}
  ")
endif()

if(CMAKE_GENERATOR MATCHES "Ninja" AND CMAKE_VERSION VERSION_GREATER_EQUAL 3.27.0 AND CMAKE_VERSION VERSION_LESS 3.27.8)
  message(WARNING "CMake 3.27.0..3.27.7 has a bug with Ninja causing build failures. Suggest using CMake outside this range or
  cmake -Bbuild -G \"Unix Makefiles\"
  ")
endif()

# --- summary
set_package_properties(GLOW PROPERTIES
TYPE OPTIONAL
URL "https://www2.hao.ucar.edu/modeling/glow/code"
DESCRIPTION "NCAR GLOW model")
#PURPOSE "Gemini uses GLOW for modeling of auroral emissions vs. wavelength.")

set_package_properties(MPI PROPERTIES
TYPE REQUIRED
DESCRIPTION "GEMINI3D MPI standard parallelization")
#PURPOSE "MPI gives massively parallel computation")

set_package_properties(MUMPS PROPERTIES
TYPE REQUIRED
URL "https://mumps-solver.org/"
DESCRIPTION "parallel direct sparse solver"
PURPOSE "MUMPS solves potential")

set_package_properties(SCALAPACK PROPERTIES
TYPE REQUIRED
URL "http://www.netlib.org/scalapack/"
DESCRIPTION "parallel linear algebra"
PURPOSE "MUMPS solves potential in parallel using Scalapack")

set_package_properties(LAPACK PROPERTIES
TYPE REQUIRED
URL "http://www.netlib.org/lapack/"
DESCRIPTION "linear algebra library"
PURPOSE "LAPACK solves parabolic and elliptical partial differential equations")

# --- options

add_feature_info(DevMode dev "Gemini developer mode")

add_feature_info(GLOW glow "airglow / aurora model")
add_feature_info(HWM14 hwm14 "HWM14 neutral winds model")

add_feature_info(PyGemini python "simulation generation, HPC script generator and plotting")
add_feature_info(MatGemini matlab "checks not as extensive as Python, and slow")

# print to screen
feature_summary(WHAT ENABLED_FEATURES DISABLED_FEATURES)
