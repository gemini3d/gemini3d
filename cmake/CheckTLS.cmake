cmake_minimum_required(VERSION 3.19)  # just FYI

function(check_tls)
# some CMake may not have SSL/TLS enabled, or may have missing/broken system certificates.
# this is a publicly-usable service (as per their TOS)

set(url https://www.howsmyssl.com/a/check)
set(temp ${PROJECT_BINARY_DIR}/test_ssl.json)

if(NOT EXISTS ${temp})
  file(DOWNLOAD ${url} ${temp} INACTIVITY_TIMEOUT 5)
endif()

file(READ ${temp} json)

string(JSON rating ERROR_VARIABLE e GET ${json} rating)

message(VERBOSE "TLS status: ${rating}")
if(NOT rating)
  message(WARNING "TLS seems to be broken on your system. Download will probably fail.  ${rating}")
endif()

endfunction(check_tls)
