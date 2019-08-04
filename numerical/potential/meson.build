# FIXME: debug Meson internals
# mumps = dependency('MUMPS', cmake_module_path : 'cmake/Modules')
mumps_root = get_option('MUMPS_ROOT')
scalapack_root = get_option('SCALAPACK_ROOT')

if mumps_root == ''
  mumpsinc = '/usr/include'
  mumpslib = fc.find_library('dmumps', required : false)
else
  message('searching for MUMPS in ' + mumps_root)
  mumpsinc = mumps_root / 'include'
  mumpslib = [fc.find_library('dmumps', required : false, dirs : mumps_root),
              fc.find_library('mumps_common', required : false, dirs : mumps_root)]
endif
if mumpslib.found()
  mumps = declare_dependency(include_directories: mumpsinc, dependencies: mumpslib)
else
  mumps_proj = subproject('mumps')
  mumps = mumps_proj.get_variable('mumps')
endif

if scalapack_root == ''
  scalapack = fc.find_library('scalapack-openmpi', required: false)
  if not scalapack.found()
    scalapack = fc.find_library('scalapack')
  endif
else
  message('searching for SCALAPACK in ' + scalapack_root)
  scalapack = fc.find_library('scalapack', dirs : scalapack_root / 'lib')
endif

blacs = fc.find_library('blacs-openmpi', required: false)

threads = dependency('threads')

potential = library('potential',
  sources: ['potential_comm_mumps.f90', 'potential_mumps.F90', 'boundary_conditions/potentialBCs_mumps.f90'],
  link_with: [calculus, collisions, fsutils, mpimod, const, grid, timeutils],
  dependencies: [mumps, scalapack, blacs, mpi, threads],
  fortran_args: real_bits)
