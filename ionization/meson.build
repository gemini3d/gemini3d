fang = library('ionize_fang', 'fang.f90',
  link_with: const)

fang_run = library('fang_run', 'fang_run.f90',
  link_with: [fang, msis, const])

test_fang_exe = executable('test_fang', 'test_fang.f90',
  link_with: [fang_run, const])

test('Fang', python,
  args: ['test_fang.py', test_fang_exe],
  timeout: 10,
  workdir: meson.current_source_dir(),
  suite: 'unit')

precipBCs = library('precipBCs', 'boundary_conditions/precipBCs_mod.f90',
  link_with: [mpimod, const, reader, grid, mesh, interp, timeutils],
  dependencies: mpi)


# ===== GLOW library
ion_libs = [const, neutral, grid, mesh, timeutils, mpimod, fang]
if get_option('glow').enabled()
  glow_proj = subproject('cglow')
  ion_libs += glow_proj.get_variable('cglow')

  ionization = library('ionization_glow',
    sources: ['ionization.f90', 'glow_run.F90'],
    link_with: ion_libs,
    dependencies : mpi,
    fortran_args: '-DDATADIR="' + glow_proj.get_variable('data_dir') + '/"')
else
  ionization = library('ionization',
    sources: ['ionization.f90', 'glow_dummy.f90'],
    link_with: ion_libs,
    dependencies : mpi)
endif
