timeutils = library('timeutils', 'timeutils.f90',
  link_with: const)

test_rollover = executable('test_rollover', 'test_dayrollover.f90',
  link_with: [timeutils, const])
test('Rollover', test_rollover)

test_sza = executable('test_sza', 'test_sza.f90',
  link_with: [timeutils, const])
test('SZA', test_sza)