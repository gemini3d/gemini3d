function model_setup_interp(p)
narginchk(1,1)
validateattributes(p, {'struct'}, {'scalar'}, mfilename, 'parameters', 1)

%% check directories
simdir_inputs = [p.simdir, '/inputs'];
eqdir_inputs = [p.eqdir, '/inputs'];

%ADD PATHS FOR FUNCTIONS
cwd = fileparts(mfilename('fullpath'));
gemdir = [cwd, '/..'];
addpath([gemdir, '/script_utils'],[gemdir, '/setup'],[gemdir, '/setup/gridgen'],[gemdir, '/vis'])

%% copy equilibrium to new directory so Fortran can read and upsample
makedir(simdir_inputs)
copyfile(eqdir_inputs, simdir_inputs)

%% GRID GENERATION
xg = makegrid_cart_3D(p);

eq2dist(p, xg);

%% potential boundary conditions
if isfield(p, 'Efield_fracwidth')
  if p.lxp == 1 || p.lyp == 1
    Efield_BCs_2d(p);
  else % 3D
    Efield_BCs_3d(p);
  end
end
%% aurora
if isfield(p, 'precip_latwidth') && isfield(p, 'precip_lonwidth')
  if p.lxp == 1 || p.lyp == 1
    particles_BCs_2d(p)
  else
    particles_BCs_3d(p)
  end
end

end % function
