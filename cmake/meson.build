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

# check if MPI compiler wrapper is usable
if fc.run('use mpi; end', name: 'Fortran MPI compiler wrapper').returncode() == 0
  mpi = []
  mpiroot = get_option('MPI_ROOT')
  if mpiroot == ''
    error('please specify -DMPI_ROOT= root path of MPI lib e.g. ~/lib/openmpi-3.1.4')
  endif
  mpiexec = find_program(mpiroot / 'bin/mpiexec')
else
  mpi = dependency('mpi', language : 'fortran')
  mpiret = fc.run('use mpi; end', dependencies : mpi, name: 'Fortran MPI')
  if mpiret.returncode() != 0
    error('MPI does not link properly: ' + mpiret.stderr() + ' ' + mpiret.stdout())
  endif
  mpiexec = find_program('mpiexec')
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


pymod = import('python')
python = pymod.find_installation('python3')
if python.language_version().version_compare('< 3.6')
  error('Python ' + python.language_version() + ' is too old, you must have at least Python 3.6. Try "conda activate" or "module load python3" or similar.')
endif

is_ci = run_command(python, meson.source_root() / 'script_utils/meson_get_environment_variable.py', 'CI').stdout().to_lower() == 'true'

python_test = python
if run_command(python, '-c', 'import gemini').returncode() != 0
  message('Need to install PyGemini by "pip install -e ' + meson.source_root() + '"')
  python_test = disabler()
else
  message('Using Python ' + python.language_version())
endif

octave = find_program('octave', required: false, disabler: true)

matlab = disabler()
if get_option('usematlab')  # buggy for Matlab < R2019a, may hang "meson setup" silently for two minutes
  matlab = find_program('matlab', required: false, disabler: true)
  if run_command(python, meson.source_root() / 'script_utils/meson_matlab_has_batch.py').returncode() != 0
    matlab = disabler()
  endif
endif