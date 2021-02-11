include(FeatureSummary)

# --- warnings

if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)

  if(WIN32)
    set(_min 3.19.2)
  else()
    set(_min 3.19.3)
  endif()

  if(CMAKE_VERSION VERSION_LESS ${_min})
    message(WARNING "Intel compilers may not work properly with CMake < ${_min}")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL IntelLLVM)
  set(_min 3.20)

  if(CMAKE_VERSION VERSION_LESS ${_min})
    message(WARNING "Intel LLVM compilers may not work properly with CMake < ${_min}")
  endif()
endif()

# --- recommendations

if(NOT CMAKE_GENERATOR MATCHES Ninja)
  message(STATUS "Ninja builds/rebuilds much faster than Make for any software project.
  Recommendation: Install Ninja build system:
      cmake -P ${PROJECT_SOURCE_DIR}/scripts/install_ninja.cmake")
endif()

if(CMAKE_VERSION VERSION_LESS 3.17)
  message(STATUS "Recommendation: Install CMake >= 3.17 by:
    cmake -P ${PROJECT_SOURCE_DIR}/scripts/install_cmake.cmake")
endif()

if(mpi AND NOT HWLOC_FOUND)
  message(STATUS "Recommendation: consider installing HWLOC for auto-detect CPU count from gemini3d.run:
    cmake -P ${PROJECT_SOURCE_DIR}/scripts/install_hwloc.cmake")
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
DESCRIPTION "NCAR GLOW model"
PURPOSE "Gemini uses GLOW for modeling of auroral emissions vs. wavelength.")

set_package_properties(HDF5 PROPERTIES
TYPE RECOMMENDED
# URL "https://www.hdfgroup.org/solutions/hdf5"
# DESCRIPTION "HDF5 file I/O"
PURPOSE "Gemini uses NASA standard format HDF5 files to read and write compressed data.")

set_package_properties(MPI PROPERTIES
TYPE RECOMMENDED
DESCRIPTION "GEMINI MPI-2 standard parallelization"
PURPOSE "MPI gives massively parallel computation")

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

set_package_properties(Python3 PROPERTIES
TYPE OPTIONAL
# URL "http://www.python.org/"
# DESCRIPTION "Python runtime"
PURPOSE "PyGemini is the standard user interface for Gemini input/output/plotting")

set_package_properties(HWLOC PROPERTIES
TYPE RECOMMENDED
URL "https://www.open-mpi.org/projects/hwloc/"
DESCRIPTION "portable abstraction of CPU hierarchical topology"
PURPOSE "Determine the number of physical CPU cores on the host computer for gemini3d.run parallel run")

# --- options

add_feature_info(DevMode dev "Gemini developer mode")
add_feature_info(MPI mpi "GEMINI MPI-2 standard parallelization")

add_feature_info(GLOW glow "airglow / aurora model")
add_feature_info(HWM14 hwm14 "HWM14 neutral winds model")
add_feature_info(MSIS2.0 msis20 "NRL MSIS 2.0 neutral atmosphere model")

add_feature_info(NetCDF4 netcdf "file read / write")

add_feature_info(PyGemini PYGEMINI_DIR "simulation generation, HPC script generator and plotting")
add_feature_info(MatGemini MATGEMINI_DIR "checks not as extensive as Python, and slow")

add_feature_info(HDF5 hdf5 "file read / write")
add_feature_info(AutoHDF5 hdf5_external "build HDF5")

add_feature_info(AutoMumps mumps_external "build Mumps")
add_feature_info(AutoScalapack scalapack_external "build Scalapack")
add_feature_info(AutoLapack lapack_external "build Lapack")

# print to screen
feature_summary(WHAT ENABLED_FEATURES)
# PACKAGES_FOUND)
