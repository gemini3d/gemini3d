fc = meson.get_compiler('fortran')

os = host_machine.system()
if os == 'windows' and fc.get_id() != 'intel'
  error('On native Windows, only the Intel compilers (ifort) works at this time, because OpenMPI is not available for Windows. Try Windows Subsystem for Linux.')
endif

# -- compiler options
quiet = fc.first_supported_argument(['-w', '-nowarn', '/nowarn'])
impnone = fc.first_supported_argument(['-fimplicit-none', '-warn declarations', '/warn:declarations', '-Mdclchk'])
std18 = fc.first_supported_argument(['-std=f2018', '-stand f18', '/stand:f18'])
bounds = fc.first_supported_argument(['-Werror=array-bounds', '-check bounds', '/check:bounds', '-Mbounds'])
if fc.get_id() == 'gcc'
  add_project_arguments('-Wno-unused-variable', '-Wno-maybe-uninitialized', '-Wno-unused-dummy-argument', '-Wno-unused-function',
    impnone, std18, language : 'fortran')
  if buildtype.contains('debug')
    add_project_arguments(bounds, '-fcheck=all', language: 'fortran')
  endif
elif fc.get_id() == 'intel'
  add_project_arguments(std18, impnone, language: 'fortran')
  if buildtype.contains('debug')
    # -heap-arrays didn't appear to be needed on Linux?
    add_project_arguments(bounds, language: 'fortran')
  endif
elif fc.get_id() == 'intel-cl'
  # /fpp is necessary for any Fortran preprocessing statement with Intel ifort on Windows
  # https://software.intel.com/en-us/fortran-compiler-developer-guide-and-reference-fpp
   add_project_arguments('/fpp', std18, impnone, language: 'fortran')
   if buildtype.contains('debug')
    add_project_arguments(bounds, '/heap-arrays', language: 'fortran')
  endif
elif ['pgi', 'flang'].contains(fc.get_id())
  add_project_arguments(impnone, language: 'fortran')
  if buildtype.contains('debug')
    add_project_arguments(bounds, language: 'fortran')
  endif
endif

# if MPI_ROOT, assume compiler wrapper is usable
mpiroot = get_option('MPI_ROOT')
if mpiroot == ''
  mpiexec = find_program('mpiexec')
  mpi = dependency('mpi', language : 'fortran')
  mpiret = fc.run('use mpi; end', dependencies : mpi, name: 'Fortran MPI')
  if mpiret.returncode() != 0
    warning('MPI_Fortran does not link properly: ' + mpiret.stderr() + ' ' + mpiret.stdout())
  endif
else
  fc.run('use mpi; integer :: i; call mpi_init(i); call mpi_finalize(i); end', name: 'Fortran MPI wrapper check')
  mpiexec = find_program(mpiroot / 'mpiexec')
  mpi = []
endif


# we don't use fallback because some systems have screwed-up libs factory installed
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

pymod = import('python')
python = pymod.find_installation('python3')
if python.language_version().version_compare('< 3.6')
  error('Python ' + python.language_version() + ' is too old, you must have at least Python 3.6. Try "conda activate" or "module load python3" or similar.')
endif

# not needed at the moment
# is_ci = run_command(python, meson.source_root() / 'script_utils/meson_get_environment_variable.py', 'CI').stdout().to_lower() == 'true'

python_test = python
if run_command(python, '-c', 'import gemini').returncode() != 0
  message('Need to install PyGemini by "pip install -e ' + meson.source_root() + '"')
  python_test = disabler()
else
  message('Using Python ' + python.language_version())
endif

octave = find_program('octave', required: false, disabler: true)
# must do .found check since run_command(disabler()) has undefined returncode.
if octave.found() and run_command(octave, '--eval', 'assert(exist("validateattributes", "file")==2)').returncode() != 0
  message('Not using too-old Octave for comparison tests.')
  octave = disabler()
endif

matlab = disabler()
if get_option('usematlab')  # buggy for Matlab < R2019a, may hang "meson setup" silently for two minutes
  matlab = find_program('matlab', required: false, disabler: true)
  # must do .found check since run_command(disabler()) has undefined returncode.
  if matlab.found() and run_command(python, meson.source_root() / 'script_utils/meson_matlab_has_batch.py').returncode() != 0
    matlab = disabler()
  endif
endif
