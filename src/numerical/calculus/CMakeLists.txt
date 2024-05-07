configure_file(calculus.in.f90 calculus.f90 @ONLY)
add_library(calculus OBJECT
${CMAKE_CURRENT_BINARY_DIR}/calculus.f90 gradient.f90 integral.f90 div.f90
)
target_link_libraries(calculus PRIVATE meshobj)
