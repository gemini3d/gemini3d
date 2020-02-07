io_src = files(
'io.f90',
'aurora.f90',
'mag.f90',
'plasma.f90',
'output.f90',
'input.f90'
)

# accomodate intel compiler non-standard inquire(directory=)
if ['intel', 'intel-cl'].contains(fc.get_id())
  io_src += files('path_exists_intel.f90')
else
  io_src += files('path_exists.f90')
endif

# FIXME: when validated, use if hdf5_interface.found() as for output.
io_src += files('plasma_input_raw.f90')

if hdf5_interface.found()
  io_src += files('plasma_output_hdf5.f90', 'mag_hdf5.f90', 'aurora_hdf5.f90')
elif netcdf_interface.found()
  io_src += files('plasma_output_ncdf.f90', 'mag_ncdf.f90', 'aurora_ncdf.f90')
else
  io_src += files('plasma_output_raw.f90', 'mag_raw.f90', 'aurora_raw.f90')
endif

# === reader

if hdf5_interface.found()
  reader_src = files('reader_hdf5.f90')
  reader_dep = hdf5_interface
elif netcdf_interface.found()
  reader_src = files('reader_nc4.f90')
  reader_dep = netcdf_interface
else
  reader_src = files('reader_raw.f90')
  reader_dep = []
endif
reader = library('reader',
  sources: ['reader.f90', reader_src],
  link_with: [const],
  dependencies: reader_dep)

# == other

pathlib_src = files('pathlib.F90')

realpath_src = files('io_dummy.f90')
if os == 'windows'
  code = '''
#include <stdlib.h>
#include <string.h>
int main(void) {
char resolved_path[_MAX_PATH];
_fullpath(resolved_path, "..", _MAX_PATH);
if(strcmp("..", resolved_path) == 0 || strcmp(".", resolved_path) == 0) return 1;
return 0; }
'''
  if (cc.run(code, name: 'Windows _fullpath').returncode() == 0)
    realpath_src = files('io_windows.f90')
  endif
elif os != 'darwin'  # FIXME: make test in Fortran not C
  code = '''
#include <limits.h>
#include <stdlib.h>
#include <string.h>
int main(void) {
char resolved_path[PATH_MAX];
realpath("..", resolved_path);
if(strcmp("..", resolved_path) == 0 || strcmp(".", resolved_path) == 0) return 1;
return 0; }
'''
  if (cc.run(code, name: 'POSIXrealpath').returncode() == 0)
    realpath_src = files('io_posix.f90')
  endif
endif


pathlib = library('pathlib',
  sources: [pathlib_src, realpath_src])

logging = library('logging', 'logging.f90')
