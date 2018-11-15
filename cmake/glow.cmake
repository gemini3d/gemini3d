function(add_glow)
  # FIXME: this will be handled far more elegantly upon merge.
  
  find_package(Git)
  if(Git_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init)
  endif()

  add_library(cglow ionization/glow/glow.f90 ionization/glow/cglow.f90 
      ionization/glow/fieldm.f ionization/glow/solzen.f90 ionization/glow/ssflux.f90
      ionization/glow/rcolum.f90 ionization/glow/qback.f90
      ionization/glow/etrans.f90 ionization/glow/exsect.f
      ionization/glow/gchem.f90 ionization/glow/bands.f90 ionization/glow/qproton.f90
      ionization/glow/ephoto.f90
      ionization/glow/egrid.f90 ionization/glow/maxt.f90)
  if(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
    target_compile_options(cglow PRIVATE -w -fno-implicit-none)
  endif()

  add_library(glow ionization/glow_run.f90)
  target_link_libraries(glow PRIVATE cglow const)

  add_library(ionization ionization/ionization.f90)
  target_link_libraries(ionization PRIVATE neutral glow)
endfunction(add_glow)

add_glow()
