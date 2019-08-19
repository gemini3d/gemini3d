fc = meson.get_compiler('fortran')
fcid = fc.get_id()
# -- compiler options
if fcid == 'gcc'
  add_project_arguments('-Wno-unused-variable', '-Wno-maybe-uninitialized', '-Wno-unused-dummy-argument', '-Wno-unused-function',
  '-fimplicit-none',
  language : 'fortran')
  if fc.version().version_compare('>=8.0')
    add_project_arguments('-std=f2018', language: 'fortran')
  endif
  if buildtype.contains('debug')
    add_project_arguments('-Werror=array-bounds', '-fcheck=all', language: 'fortran')
  endif
elif fcid == 'intel'
  add_project_arguments('-stand', 'f18', '-warn declarations', language: 'fortran')
  if buildtype.contains('debug')
    # -heap-arrays didn't appear to be needed on Linux?
    add_project_arguments('-check bounds', language: 'fortran')
  endif
elif fcid == 'intel-cl'
  # /fpp is necessary for any Fortran preprocessing statement with Intel ifort on Windows
  # https://software.intel.com/en-us/fortran-compiler-developer-guide-and-reference-fpp
   add_project_arguments('/fpp', '/stand:f18', '/warn:declarations', language: 'fortran')
   if buildtype.contains('debug')
    add_project_arguments('/check:bounds', '/heap-arrays', language: 'fortran')
  endif
elif fcid == 'pgi'
  add_project_arguments('-Mdclchk', language: 'fortran')
  if buildtype.contains('debug')
    add_project_arguments('-Mbounds', language: 'fortran')
  endif
endif
# need find_library to search specific dirs at this time
mpi_root = get_option('MPI_ROOT')
if mpi_root == ''
  mpi = dependency('mpi', language : 'fortran', required : true)
else
  message('searching for MPI in ' + mpi_root)
  mpi = fc.find_library('mpi', required : true, dirs : mpi_root / 'lib')
endif

mpiret = fc.run('use mpi; end', dependencies : mpi, name: 'Fortran MPI')
if mpiret.returncode() != 0
  error('MPI does not link properly: ' + mpiret.stderr() + mpiret.stdout() + mpi_root + ' ' + mpi.name() + ' ' + mpi.version())
endif

system_blas = get_option('system_blas')
lapack_root = get_option('LAPACK_ROOT')
if system_blas
  if lapack_root == ''
    lapack = dependency('lapack', required: false)
    blas = dependency('blas', required: false)
    if not blas.found()  # necessary for some systems incl. CentOS 7
      blas = fc.find_library('blas', required: false)
    endif
  else
    message('searching for LAPACK & BLAS in ' + lapack_root)
    lapack = fc.find_library('lapack', dirs : lapack_root, required: false)
    blas = fc.find_library('blas', dirs : lapack_root, required: false)
  endif
else
  lapack = disabler()
  blas = disabler()
endif
if not lapack.found() or not blas.found()
  lapack_proj = subproject('lapack')
  lapack = lapack_proj.get_variable('lapack')
  blas = lapack_proj.get_variable('blas')
endif

if use_hdf  # see HDF5 branch, not yet fully ready.
 hdf5 = dependency('hdf5')
endif