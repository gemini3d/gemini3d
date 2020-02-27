function model_setup_interp(p)
narginchk(1,1)
validateattributes(p, {'struct'}, {'scalar'}, mfilename, 'parameters', 1)

%% ADD PATHS FOR FUNCTIONS
cwd = fileparts(mfilename('fullpath'));
gemdir = [cwd, '/..'];
addpath([gemdir, '/script_utils'],[gemdir, '/setup/gridgen'],[gemdir, '/vis'])

%% GRID GENERATION
xg = makegrid_cart_3D(p);

eq2dist(p, xg);

%% potential boundary conditions
if isfield(p, 'flagE0file') && p.flagE0file
  if p.lxp == 1 || p.lyp == 1
    Efield_BCs_2d(p);
  else % 3D
    Efield_BCs_3d(p);
  end
end

%% aurora
if isfield(p, 'flagprecfile') && p.flagprecfile
  particles_BCs(p)
end

end % function
