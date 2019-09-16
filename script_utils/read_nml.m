function params = read_nml(filename)

narginchk(1,1)

params = read_nml_group(filename, 'base');
if params.flagdneu
  params = merge_struct(params, read_nml_group(filename, 'neutral_perturb'));
end
if params.flagprecfile
  params = merge_struct(params, read_nml_group(filename, 'precip'));
end
if params.flagE0file
  params = merge_struct(params, read_nml_group(filename, 'efield'));
end
if params.flagglow
  params = merge_struct(params, read_nml_group(filename, 'glow'));
end

end % function