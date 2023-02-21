function(patch_lapack95)

if(realbits EQUAL 32)
  set(sgbsv "call sgbsv(N,LKL,KU,NRHS,A,LDA,LPIV,B,N,LINFO)")
  set(dgbsv "B = 0
  error stop 'LAPACK95:gbsv Gemini3D is in real32 mode. File GitHub issue for Gemini'")
else()
  set(dgbsv "call dgbsv(N,LKL,KU,NRHS,A,LDA,LPIV,B,N,LINFO)")
  set(sgbsv "B = 0
  error stop 'LAPACK95:gbsv Gemini3D is in real64 mode. File GitHub issue for Gemini'")
endif()

configure_file(gbsv.in.f90 gbsv.f90 @ONLY)

endfunction(patch_lapack95)

patch_lapack95()

add_library(gbsv OBJECT ${CMAKE_CURRENT_BINARY_DIR}/gbsv.f90)
