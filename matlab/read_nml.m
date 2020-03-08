function p = read_nml(path)

% for reading simulation config*.nml. Fortran namelist is a standard
% format.

narginchk(1,1)

filename = get_configfile(path);

p = read_nml_group(filename, 'base');
p.indat_file = absolute_path(p.indat_file);
p.indat_size = absolute_path(p.indat_size);
p.indat_grid = absolute_path(p.indat_grid);
if ~isfield(p, 'simdir')
  p.simdir = fileparts(p.indat_size);
end

if ~isfield(p, 'nml')
  p.nml = filename;
end

try %#ok<TRYNC>
p = merge_struct(p, read_nml_group(filename, 'setup'));
end
if isfield(p, 'eqdir')
  p.eqdir = absolute_path(p.eqdir);
end

if isfield(p, 'flagdneu') && p.flagdneu
  p = merge_struct(p, read_nml_group(filename, 'neutral_perturb'));
end
if ~isfield(p, 'mloc')
  p.mloc=[];
end

if isfield(p, 'flagprecfile') && p.flagprecfile
  p = merge_struct(p, read_nml_group(filename, 'precip'));
  p.prec_dir = absolute_path(p.prec_dir);
end

if isfield(p, 'flagE0file') && p.flagE0file
  p = merge_struct(p, read_nml_group(filename, 'efield'));
  p.E0_dir = absolute_path(p.E0_dir);
end

if isfield(p, 'flagglow') && p.flagglow
  p = merge_struct(p, read_nml_group(filename, 'glow'));
end

end % function
