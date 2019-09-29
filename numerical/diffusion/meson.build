parabolic = library('PDEparabolic',
  sources: ['PDEparabolic.f90', '../../vendor/lapack95/gbsv.f90'],
  link_with: const,
  dependencies: lapack)

diffusion = library('diffusion', 'diffusion.f90',
  link_with: [const, grid, parabolic])

test_diffusion1D_exe = executable('test_diffusion1D', 'test_diffusion1D.f90',
  link_with: [const, diffusion, parabolic])
test('diffusion1D', test_diffusion1D_exe,
  timeout: 5,
  suite: 'unit')

test('Diffusion1D_octave', octave,
  args: ['-q', '--eval', 'test_diffusion1D("' + meson.build_root() / 'test_diffusion1d.dat")'],
  timeout: 10,
  workdir: meson.current_source_dir(),
  suite: 'unit')

# must have .current_source_dir() for Python script
test('Diffusion1D_python', python,
  args: [meson.current_source_dir() / 'test_diffusion1D.py', meson.build_root() / 'test_diffusion1d.dat'],
  timeout: 10,
  suite: 'unit')