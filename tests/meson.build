mpiexec = find_program('mpiexec')
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

test_formats = executable('test_formats', 'test_formats.f90',
  link_with: [fsutils, timeutils])
test('Formats', test_formats, suite : 'unit')


octave = find_program('octave', required: false)

matlab = disabler()
if get_option('usematlab')  # buggy for Matlab < R2019a, may hang "meson setup" silently for two minutes
  matlab = find_program('matlab', required: false)
  if run_command(python, meson.source_root() / 'script_utils/meson_matlab_has_batch.py').returncode() != 0
    matlab = disabler()
  endif
endif