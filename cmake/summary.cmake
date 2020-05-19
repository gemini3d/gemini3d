include(FeatureSummary)

set_package_properties(MPI PROPERTIES
    URL "https://www.open-mpi.org/"
    DESCRIPTION "OpenMPI, IntelMPI, MPICH and MS-MPI are known to work with GEMINI"
    PURPOSE "MPI is essential to GEMINI for massively parallel computation.")

add_feature_info(WithoutMPI no_mpi "partial build of unit tests only")

add_feature_info(GLOW glow "airglow / aurora model")
add_feature_info(HDF5 HDF5OK "file read / write")
add_feature_info(NetCDF4 NCDFOK "file read / write")
add_feature_info(Python python_ok "simulation generation, HPC script generator and plotting")

add_feature_info(AutoMumps mumps_external "auto-build Mumps")
add_feature_info(AutoScalapack scalapack_external "auto-build Scalapack")
add_feature_info(AutoLapack lapack_external "auto-build Lapack")

# print to screen
feature_summary(WHAT ENABLED_FEATURES DISABLED_FEATURES)