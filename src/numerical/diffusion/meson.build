parabolic = library('PDEparabolic',
  sources: ['PDEparabolic.f90', '../../vendor/lapack95/gbsv.f90'],
  link_with: const,
  dependencies: lapack)

test_diffusion1D_exe = executable('test_diffusion1D', 'test_diffusion1D.f90',
  link_with: [const, parabolic],
  dependencies: hdf5_interface)
test('diffusion1D', test_diffusion1D_exe,
  timeout: 5,
  suite: 'unit')

test('Diffusion1D_octave', octave,
  args: ['-q', '--eval', 'test_diffusion1D("' + meson.build_root() / 'test_diffusion1d.h5")'],
  timeout: 10,
  workdir: meson.current_source_dir(),
  suite: 'unit', is_parallel: false, priority: -100)

test('Diffusion1D_python', python,
  args: ['test_diffusion1D.py', meson.build_root() / 'test_diffusion1d.h5'],
  timeout: 10,
  workdir: meson.current_source_dir(),
  suite: 'unit', is_parallel: false, priority: -100)


diffusion = library('diffusion', 'diffusion.f90',
  link_with: [const, grid, mesh, parabolic])
