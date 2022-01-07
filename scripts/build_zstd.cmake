cmake_minimum_required(VERSION 3.20...3.22)

file(READ ${CMAKE_CURRENT_LIST_DIR}/versions.json _j)
string(JSON zstd_version GET ${_j} zstd)


set(prefix "~/zstd-${zstd_version}")

set(CMAKE_TLS_VERIFY true)

if(DEFINED ENV{TMPDIR})
  set(tmpdir $ENV{TMPDIR})
elseif(IS_DIRECTORY /var/tmp)
  set(tmpdir /var/tmp)
elseif(IS_DIRECTORY /tmp)
  set(tmpdir /tmp)
else()
  set(tmpdir ~/tmp)
endif()

if(CMAKE_VERSION VERSION_LESS 3.21)
  get_filename_component(tmpdir ${tmpdir} ABSOLUTE)
else()
  file(REAL_PATH ${tmpdir} tmpdir EXPAND_TILDE)
endif()

set(name zstd-${zstd_version}.tar.gz)
set(archive ${tmpdir}/${name})

set(src ${tmpdir}/zstd-${zstd_version}/build/cmake)
set(build ${src}/build)

if(NOT IS_DIRECTORY ${src})
  file(DOWNLOAD https://github.com/facebook/zstd/releases/download/v${zstd_version}/${name} ${archive}
    INACTIVITY_TIMEOUT 15)
  file(ARCHIVE_EXTRACT INPUT ${archive} DESTINATION ${tmpdir})
endif()

execute_process(COMMAND ${CMAKE_COMMAND} --install-prefix=${prefix} -S ${src} -B ${build}
COMMAND_ERROR_IS_FATAL ANY)

execute_process(COMMAND ${CMAKE_COMMAND} --build ${build} --parallel
COMMAND_ERROR_IS_FATAL ANY)

execute_process(COMMAND ${CMAKE_COMMAND} --install ${build}
COMMAND_ERROR_IS_FATAL ANY)

message(STATUS "Please add ${CMAKE_INSTALL_PREFIX}/bin to environment variable PATH")
