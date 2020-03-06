interp = library('interp',
  sources : ['interpolation.f90','interp2d.f90'],
  link_with: const)

test_interp1_exe = executable('test_interp1', 'testinterp1.f90',
  link_with: [const, interp])
test('Interp1d', test_interp1_exe,
  suite: 'unit', timeout: 15, workdir: meson.current_build_dir())

if not hdf5_interface.found()
  subdir_done()
endif

test_interp2_exe = executable('test_interp2', 'testinterp2.f90',
  link_with: [const, interp],
  dependencies: hdf5_interface)
test('Interp2d', test_interp2_exe,
  suite: 'unit', timeout: 15, workdir: meson.current_build_dir())
test('Interp2d_python', python,
  args: [files('testinterp.py'),
         meson.current_build_dir() / 'output2d.h5'],
  timeout: 30, suite: 'unit',
  is_parallel: false, priority: -100)
test('Interp2d_octave', octave,
  args: ['-q', '--eval',
         'testinterp("' + meson.current_build_dir() / 'output2d.h5")'],
  timeout: 30, workdir: meson.current_source_dir(),  suite: 'unit',
  is_parallel: false, priority: -100)

test_interp3_exe = executable('test_interp3', 'testinterp3.f90',
  link_with: [const, interp],
  dependencies: hdf5_interface)
test('Interp3d', test_interp3_exe,
  suite: 'unit', timeout: 30, workdir: meson.current_build_dir())
test('Interp3d_python', python,
  args: [meson.current_source_dir() / 'testinterp3.py',
         meson.current_build_dir() / 'output3d.h5'],
  timeout: 90, suite: 'unit', # timeout 90 for slow cpu
  is_parallel: false, priority: -100)
test('Interp3d_octave', octave,
  args: ['-q', '--eval',
         'testinterp3("' + meson.current_build_dir() / 'output3d.h5")'],
  timeout: 30, workdir: meson.current_source_dir(), suite: 'unit',
  is_parallel: false, priority: -100)
