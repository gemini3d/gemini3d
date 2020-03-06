mumpsifc_src = bits == '32' ? files('mumps_real32.f90') : files('mumps_real64.f90')

pde_elliptic = library('pde_elliptic',
  sources: ['PDEelliptic.f90', 'elliptic2d.f90', 'elliptic3d.f90', mumpsifc_src],
  dependencies: [mpi, mumps],
  link_with: const)

pot_src = files('potential_comm_mumps.f90', 'potential_root.f90', 'potential_worker.f90',
'potential_mumps.f90', 'potential2d.f90',
'boundary_conditions/potentialBCs_mumps.f90')
pot_deps = [mumps, scalapack, blacs, lapack, mpi, threads]

potential = library('potential',
  sources: pot_src,
  link_with: [pde_elliptic, calculus, collisions, const, reader,
    grid, mesh, interp, mpimod, timeutils],
  dependencies: pot_deps)

if not hdf5_interface.found()
  subdir_done()
endif

test_potential2d = executable('test_potential2d',
  sources: 'test_potential2D.f90',
  link_with: [pde_elliptic, const],
  dependencies: [blacs, mpi, hdf5_interface])
test('potential2d', test_potential2d,
  timeout: 30, env: nomalloc_env, suite: 'unit',
  workdir: meson.current_build_dir())

test('potential2d_octave', octave,
  args: ['-q', '--eval', 'test_potential2D("' + meson.current_build_dir() / 'test_potential2d.h5")'],
  timeout: 15,
  env: nomalloc_env,
  workdir: meson.current_source_dir(),
  suite: 'unit', is_parallel: false, priority: -100)

test('potential2d_python', python,
  args: ['test_potential2d.py', meson.current_build_dir() / 'test_potential2d.h5'],
  timeout: 15,
  env: nomalloc_env,
  workdir: meson.current_source_dir(),
  suite: 'unit', is_parallel: false, priority: -100)

#test_potential3d = executable('test_potential3d', 'test_potential3D.f90',
#  dependencies: [mpi, mumps])
#test('potential3d', test_potential3d, timeout: 60, suite: 'unit')
