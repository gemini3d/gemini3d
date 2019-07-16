ionization_src = files('ionization.f90', 'boundary_conditions/precipBCs_mod.f90')
if GLOW
  ionization_src += 'glow_run.f90'
  glow_proj = subproject('cglow')
  cglow = glow_proj.get_variable('cglow')
  ionization_deps = cglow
else
  ionization_src += 'glow_dummy.f90'
  ionization_deps = []
endif

ionization = library('ionization', ionization_src,
  link_with: [ionization_deps, io, const, calculus, grid, neutral, timeutils, mpimod],
  dependencies : mpi)