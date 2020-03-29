include(FeatureSummary)

set_package_properties(MPI PROPERTIES
    URL "https://www.open-mpi.org/"
    DESCRIPTION "OpenMPI, IntelMPI, MPICH and MS-MPI are known to work with GEMINI"
    PURPOSE "MPI is essential to GEMINI for massively parallel computation.")

add_feature_info(GLOW glow "airglow / aurora model")
add_feature_info(HDF5 hdf5 "file read / write")
add_feature_info(NetCDF4 netcdf "file read / write")
add_feature_info(Python python_ok "simulation generation, HPC script generator and plotting")

# print to screen
feature_summary(WHAT ENABLED_FEATURES DISABLED_FEATURES)