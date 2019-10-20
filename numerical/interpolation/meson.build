interp = library('interp', 'interpolation.f90',
  link_with: const)

test_interp1_exe = executable('test_interp1', 'testinterp1.f90',
  link_with: [const, interp])
test('Interp1d', test_interp1_exe,
  suite: 'unit', timeout: 15)

test_interp2_exe = executable('test_interp2', 'testinterp2.f90',
  link_with: [const, interp])
test('Interp2d', test_interp2_exe,
  suite: 'unit', timeout: 15)

test_interp3_exe = executable('test_interp3', 'testinterp3.f90',
  link_with: [const, interp])


#== Octave
test('Interp2d_octave', octave,
  args: ['-q', '--eval',
         'testinterp("' + meson.build_root() / 'output2D.dat",' + get_option('realbits') +')'],
  timeout: 30,
  workdir: meson.current_source_dir(),
  suite: 'unit')


#== Python
# must have .current_source_dir() for Python script
test('Interp2d_python', python,
  args: [meson.current_source_dir() / 'testinterp.py',
         meson.build_root() / 'output2D.dat',
         get_option('realbits')],
  timeout: 30,
  suite: 'unit')

# interp3d not parallel because of 600+ MB RAM usage (for small RAM systems)
test('Interp3d_python', python,
  args: [meson.current_source_dir() / 'testinterp3.py',
         meson.build_root() / 'output3D.dat',
         get_option('realbits'),
         test_interp3_exe],
  timeout: 90, # need 90 for slow CPU
  suite: 'unit',
  is_parallel: false)
