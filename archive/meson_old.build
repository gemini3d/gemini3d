octave = find_program('octave', required : false)
matlab = find_program('matlab', required : false)
if matlab.found()
  octmat = matlab
  evalbatch = '-batch'
elif octave.found()
  octmat = octave
  evalbatch = '--eval'
else
  octmat = []
  evalbatch = ''
  message('Neither Matlab or Octave were found, "Compare" tests will fail.')
endif



  if octmat.found()
    octmat_args = [evalbatch, 'compare_all("' + test_dir + '","' + zenodo[i]['dir'] + '")']
    test('Compare' + i, octmat, args : octmat_args,
      suite : i, timeout: 60, workdir : meson.source_root() / 'tests')
  endif