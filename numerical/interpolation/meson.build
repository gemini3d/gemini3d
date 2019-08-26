interp = library('interp', 'interpolation.f90',
  link_with: const)

test_interp1_exe = executable('test_interp1', 'testinterp1.f90',
  link_with: [const, interp])
test('Interp1d', test_interp1_exe,
  suite: 'unit')

test_interp2_exe = executable('test_interp2', 'testinterp2.f90',
  link_with: [const, interp])
test('Interp2d', test_interp2_exe,
  suite: 'unit')

test_interp3_exe = executable('test_interp3', 'testinterp3.f90',
  link_with: [const, interp])


#== Octave
test('Interp2d_octave', octave,
  args: ['-q', '--eval', 'testinterp("' + meson.build_root() / 'output2D.dat")'],
  timeout: 10,
  workdir: meson.current_source_dir(),
  suite: 'unit')


#== Python
# must have .current_source_dir() for Python script
test('Interp2d_python', python,
  args: [meson.current_source_dir() / 'testinterp.py', meson.build_root() / 'output2D.dat'],
  timeout: 10,
  suite: 'unit')

test('Interp3d_python', python,
  args: [meson.current_source_dir() / 'testinterp3.py', meson.build_root() / 'output3D.dat', test_interp3_exe],
  timeout: 10,
  suite: 'unit')
