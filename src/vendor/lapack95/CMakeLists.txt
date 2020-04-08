if(realbits EQUAL 32)
  set(sgbsv "call sgbsv(N,LKL,KU,NRHS,A,LDA,LPIV,B,N,LINFO)")
else()
  set(sgbsv "error stop 'Lapack did not seem to have real32 enabled. File GitHub issue for Gemini'")
endif()

configure_file(gbsv.in.f90 gbsv.f90)

add_library(gbsv OBJECT ${CMAKE_CURRENT_BINARY_DIR}/gbsv.f90)
