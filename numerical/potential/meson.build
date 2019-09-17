mumps_root = get_option('MUMPS_ROOT')
if bits == '64'
  mumps_name = 'dmumps'
elif bits == '32'
  mumps_name = 'smumps'
else
  error('unknown precision ' + bits)
endif
if system_blas
  # mumps = dependency('MUMPS', cmake_module_path : 'cmake/Modules')
  if mumps_root == ''
    mumpsinc = '/usr/include'
    mumpslib = fc.find_library(mumps_name, required : false)
  else
    message('searching for MUMPS in ' + mumps_root)
    mumpsinc = mumps_root / 'include'
    mumpslib = [fc.find_library(mumps_name, required : false, dirs : mumps_root),
                fc.find_library('mumps_common', required : false, dirs : mumps_root)]
  endif
else
  mumpslib = disabler()
endif
if mumpslib.found()
  mumps = declare_dependency(include_directories: mumpsinc, dependencies: mumpslib)
else
  mumps_proj = subproject('mumps')
  mumps = mumps_proj.get_variable('mumps')
endif


scalapack_root = get_option('SCALAPACK_ROOT')
if system_blas
  # dependency('scalapack') not yet working in Meson--cmake->meson bug within Meson, yielding long meson-log.txt
  # scalapack = dependency('scalapack', cmake_module_path : 'cmake/Modules')
  if scalapack_root == ''  # system
    scalapack = fc.find_library('scalapack-openmpi', required: false)
    if not scalapack.found()
      scalapack = fc.find_library('scalapack', required: false)
    endif
  else
    message('searching for SCALAPACK in ' + scalapack_root)
    scalapack = fc.find_library('scalapack', dirs : scalapack_root / 'lib')
  endif
else
  scalapack = disabler()
endif
if not scalapack.found()
  scalapack_proj = subproject('scalapack')
  scalapack = declare_dependency(link_with: scalapack_proj.get_variable('scalapack_lib'))
endif

blacs = fc.find_library('blacs-openmpi', required: false)

threads = dependency('threads')

potential = library('potential',
  sources: ['potential_comm_mumps.f90', 'potential_mumps.F90', 'boundary_conditions/potentialBCs_mumps.f90'],
  link_with: [calculus, collisions, const, grid, interp, mpimod, timeutils],
  dependencies: [mumps, scalapack, blacs, mpi, threads],
  fortran_args: real_bits)
