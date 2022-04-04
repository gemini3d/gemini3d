file(DOWNLOAD https://connectivitycheck.gstatic.com/generate_204
TIMEOUT 5
STATUS ret
TLS_VERIFY ON
)
list(GET ret 0 status)
if(NOT status EQUAL 0)
  list(GET ret 1 msg)
  message(FATAL_ERROR "Internet connectivity test failed.
  Return code: ${status}
  Error: ${msg}"
  )
endif()
