timeutils = library('timeutils', 'timeutils.f90',
  link_with: const)

# -- optional tests

test_formats = executable('test_formats', 'test_formats.f90',
  link_with: [const, timeutils])
test('Formats', test_formats, suite : 'unit')


test_rollover = executable('test_rollover', 'test_dayrollover.f90',
  link_with: [timeutils, const])
test('Rollover', test_rollover, suite : 'unit')

test_sza = executable('test_sza', 'test_sza.f90',
  link_with: [timeutils, const])
test('SZA', test_sza, suite : 'unit')