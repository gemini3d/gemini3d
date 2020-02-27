function model_setup(p)
%% determines what kind of setup is needed and does it.

narginchk(1,1)

%% ADD PATHS FOR FUNCTIONS
cwd = fileparts(mfilename('fullpath'));
gemdir = [cwd, '/..'];
addpath([gemdir, '/script_utils'],[gemdir, '/setup/gridgen'],[gemdir, '/vis'])

%% parse input
if isa(p, 'struct')
  validateattributes(p, {'struct'}, {'scalar'}, mfilename, 'parameters', 1)
elseif isa(p, 'char')
  % path to config.nml
  validateattributes(p, {'char'}, {'vector'}, mfilename, 'parameters', 1)
  p = read_nml(p);
end


%% is this equilibrium or interpolated simulation
if isfield(p, 'eqdir')
  model_setup_interp(p)
else
  model_setup_equilibrium(p)
end

end % function