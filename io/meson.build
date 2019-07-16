fsutils = library('fsutils', 'mkdir.F90', 'expanduser.f90', 'formats.f90')

logging = library('logging', 'logging.f90')

test_formats = executable('test_formats', 'test_formats.f90',
  link_with: [fsutils, timeutils])
test('Formats', test_formats)