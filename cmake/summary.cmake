include(FeatureSummary)

set_package_properties(Git PROPERTIES
TYPE REQUIRED
URL "https://git-scm.com"
DESCRIPTION "version control"
PURPOSE "Gemini retrieves the PyGemini interface using Git, if PyGemini is not already installed.")

set_package_properties(GLOW PROPERTIES
TYPE OPTIONAL
URL "https://www2.hao.ucar.edu/modeling/glow/code"
DESCRIPTION "NCAR GLOW model"
PURPOSE "Gemini uses GLOW for modeling of auroral emissions vs. wavelength.")

set_package_properties(HDF5 PROPERTIES
TYPE RECOMMENDED
# URL "https://www.hdfgroup.org/solutions/hdf5"
# DESCRIPTION "HDF5 file I/O"
PURPOSE "Gemini uses NASA standard format HDF5 files to read and write compressed data.")

set_package_properties(MPI PROPERTIES
TYPE RECOMMENDED
DESCRIPTION "GEMINI uses MPI-2 standard"
PURPOSE "MPI gives massively parallel computation")

set_package_properties(MUMPS PROPERTIES
TYPE OPTIONAL
URL "https://mumps-solver.org/"
DESCRIPTION "parallel direct sparse solver"
PURPOSE "MUMPS solves potential")

set_package_properties(SCALAPACK PROPERTIES
TYPE OPTIONAL
URL "http://www.netlib.org/scalapack/"
DESCRIPTION "parallel linear algebra"
PURPOSE "MUMPS solves potential in parallel using Scalapack")

set_package_properties(LAPACK PROPERTIES
TYPE REQUIRED
URL "http://www.netlib.org/lapack/"
DESCRIPTION "linear algebra library"
PURPOSE "LAPACK solves parabolic and elliptical partial differential equations")

set_package_properties(Python3 PROPERTIES
TYPE RECOMMENDED
# URL "http://www.python.org/"
# DESCRIPTION "Python runtime"
PURPOSE "PyGemini is the standard user interface for Gemini input/output/plotting")

# --- options

add_feature_info(MPI_gemini mpi "Use MPI parallelization")

add_feature_info(GLOW glow "airglow / aurora model")
add_feature_info(HWM14 hwm14 "HWM14 neutral winds model")
add_feature_info(MSIS2.0 msis20 "NRL MSIS 2.0 neutral atmosphere model")

add_feature_info(NetCDF4 netcdf "file read / write")

add_feature_info(Python python_ok "simulation generation, HPC script generator and plotting")
add_feature_info(Matlab Matlab_FOUND "checks not as extensive as Python, and slow")

add_feature_info(HDF5 hdf5 "file read / write")
add_feature_info(AutoHDF5 hdf5_external "auto-build HDF5")

add_feature_info(AutoMumps mumps_external "auto-build Mumps")
add_feature_info(AutoScalapack scalapack_external "auto-build Scalapack")
add_feature_info(AutoLapack lapack_external "auto-build Lapack")

# print to screen
feature_summary(WHAT ENABLED_FEATURES PACKAGES_FOUND)
