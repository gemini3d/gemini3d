function(get_url_name name out_var)
# a priori test_name strips trailing _cpp
if(name MATCHES "_cpp$")
  string(LENGTH ${name} L)
  math(EXPR M "${L}-4")
  string(SUBSTRING ${name} 0 ${M} url_name)
else()
  set(url_name ${name})
endif()

set(${out_var} ${url_name} PARENT_SCOPE)

endfunction()
