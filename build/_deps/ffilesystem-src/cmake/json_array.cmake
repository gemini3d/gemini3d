function(json_array json key var)

string(JSON L LENGTH ${json} ${key})

math(EXPR L "${L} - 1")
set(v)
foreach(i RANGE ${L})
  string(JSON s GET ${json} ${key} ${i})
  list(APPEND v ${s})
endforeach()

set(${var} ${v} PARENT_SCOPE)

endfunction()
