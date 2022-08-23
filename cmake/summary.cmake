include(FeatureSummary)

# --- recommendations

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 7.5.0)
  message(WARNING "GCC older than 7.5.0 has bugs that are likely to cause Gemini3D (and other modern programs) to fail to build.")
endif()

if(CMAKE_GENERATOR MATCHES "Visual Studio")
  message(WARNING "Visual Studio generator ${CMAKE_GENERATOR} is not supported. Please use MinGW Makefiles or Ninja:

  cmake -G Ninja -B ${PROJECT_BINARY_DIR}

  ")
endif()

# --- summary
set_package_properties(Git PROPERTIES
TYPE REQUIRED
URL "https://git-scm.com"
DESCRIPTION "version control"
PURPOSE "Git is used to auto-download the packages comprising Gemini.")

set_package_properties(GLOW PROPERTIES
TYPE OPTIONAL
URL "https://www2.hao.ucar.edu/modeling/glow/code"
DESCRIPTION "NCAR GLOW model")
#PURPOSE "Gemini uses GLOW for modeling of auroral emissions vs. wavelength.")

set_package_properties(MPI PROPERTIES
TYPE RECOMMENDED
DESCRIPTION "GEMINI MPI-2 standard parallelization")
#PURPOSE "MPI gives massively parallel computation")

set_package_properties(MUMPS PROPERTIES
TYPE RECOMMENDED
URL "https://mumps-solver.org/"
DESCRIPTION "parallel direct sparse solver"
PURPOSE "MUMPS solves potential")

set_package_properties(SCALAPACK PROPERTIES
TYPE OPTIONAL
URL "http://www.netlib.org/scalapack/"
DESCRIPTION "parallel linear algebra"
PURPOSE "MUMPS solves potential in parallel using Scalapack")

set_package_properties(LAPACK PROPERTIES
TYPE RECOMMENDED
URL "http://www.netlib.org/lapack/"
DESCRIPTION "linear algebra library"
PURPOSE "LAPACK solves parabolic and elliptical partial differential equations")

# --- options

add_feature_info(DevMode dev "Gemini developer mode")

add_feature_info(GLOW glow "airglow / aurora model")
add_feature_info(HWM14 hwm14 "HWM14 neutral winds model")
add_feature_info(MSIS2.0 msis2 "NRL MSIS 2.x neutral atmosphere model")

add_feature_info(PyGemini python "simulation generation, HPC script generator and plotting")
add_feature_info(MatGemini matlab "checks not as extensive as Python, and slow")

# print to screen
feature_summary(WHAT ENABLED_FEATURES)
