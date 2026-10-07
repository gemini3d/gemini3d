# Build the reviewed HDF5 2.2.0 profile from the publisher's SHA256-verified source.
# cmake -Dqualification_work=/tmp/hdf5-work -Dqualification_prefix=/opt/hdf5-2.2.0
#       [-Dsource_archive=/cache/hdf5-2.2.0.tar.gz] -P scripts/build_qualification_hdf5.cmake
cmake_minimum_required(VERSION 3.25)
foreach(key IN ITEMS qualification_work qualification_prefix)
  if(NOT DEFINED ${key} OR NOT IS_ABSOLUTE "${${key}}")
    message(FATAL_ERROR "${key} must be an explicit absolute path")
  endif()
endforeach()
set(archive_sha256 1a1ab8209b35586fbc1aa279ba76d102130b95badcb20ca329587219112d8c16)
file(MAKE_DIRECTORY "${qualification_work}")
if(NOT DEFINED source_archive)
  set(source_archive "${qualification_work}/hdf5-2.2.0.tar.gz")
  if(NOT EXISTS "${source_archive}")
    file(DOWNLOAD "https://support.hdfgroup.org/releases/hdf5/2.2.0/downloads/hdf5-2.2.0.tar.gz"
      "${source_archive}" EXPECTED_HASH SHA256=${archive_sha256} TLS_VERIFY ON STATUS download_status)
    list(GET download_status 0 status_code)
    if(NOT status_code EQUAL 0)
      message(FATAL_ERROR "HDF5 source download failed: ${download_status}")
    endif()
  endif()
endif()
file(SHA256 "${source_archive}" actual_sha256)
if(NOT actual_sha256 STREQUAL archive_sha256)
  message(FATAL_ERROR "HDF5 source hash mismatch; refusing build")
endif()
set(source_dir "${qualification_work}/source-${archive_sha256}")
if(EXISTS "${source_dir}" OR EXISTS "${qualification_work}/build")
  message(FATAL_ERROR "Use a fresh qualification_work directory; existing sources/builds are not trusted")
endif()
file(MAKE_DIRECTORY "${source_dir}")
# Extract only into a fresh tree so stale or edited files cannot enter the build.
file(ARCHIVE_EXTRACT INPUT "${source_archive}" DESTINATION "${source_dir}")
execute_process(COMMAND "${CMAKE_COMMAND}" -S "${source_dir}/hdf5-2.2.0" -B "${qualification_work}/build"
  -DCMAKE_BUILD_TYPE=Release "-DCMAKE_INSTALL_PREFIX=${qualification_prefix}"
  -DBUILD_TESTING=OFF -DBUILD_STATIC_LIBS=OFF -DBUILD_SHARED_LIBS=ON
  -DHDF5_BUILD_FORTRAN=ON -DHDF5_BUILD_HL_LIB=ON -DHDF5_BUILD_CPP_LIB=OFF
  -DHDF5_BUILD_TOOLS=OFF -DHDF5_BUILD_EXAMPLES=OFF
  -DHDF5_ENABLE_ZLIB_SUPPORT=ON -DHDF5_ENABLE_SZIP_SUPPORT=OFF
  COMMAND_ERROR_IS_FATAL ANY)
execute_process(COMMAND "${CMAKE_COMMAND}" --build "${qualification_work}/build" --parallel 2
  COMMAND_ERROR_IS_FATAL ANY)
execute_process(COMMAND "${CMAKE_COMMAND}" --install "${qualification_work}/build" COMMAND_ERROR_IS_FATAL ANY)
