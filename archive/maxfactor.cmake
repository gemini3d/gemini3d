# find the maximum factor of x <= y
# result variable MAXFACTOR

function(maxfactor x y)
if(x LESS 1 OR y LESS 1)
  message(FATAL_ERROR "Factors must be positive integers")
endif()

if(y LESS x)
  set(N ${y})
else()
  set(N ${x})
endif()

foreach(i RANGE 1 ${N})
  math(EXPR Mx "${x} % ${i}")
  if(Mx EQUAL 0)
    set(MAXFACTOR ${i} PARENT_SCOPE)
  endif()
endforeach()

endfunction(maxfactor)
