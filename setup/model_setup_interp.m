function [state, E] = model_setup_interp(p)
narginchk(1,1)
validateattributes(p, {'struct'}, {'scalar'}, mfilename, 'parameters', 1)

%% check directories
simdir_inputs = [p.simdir, '/inputs'];
eqdir_inputs = [p.eqdir, '/inputs'];

%ADD PATHS FOR FUNCTIONS
cwd = fileparts(mfilename('fullpath'));
gemdir = [cwd, '/../../gemini'];
assert(isfolder(gemdir), [gemdir, ' not found'])
for d = {'script_utils', 'setup', 'setup/gridgen', 'vis'}
  addpath([gemdir, filesep, d{:}])
end

%% copy equilibrium to new directory so Fortran can read and upsample
makedir(simdir_inputs)
copyfile(eqdir_inputs, simdir_inputs)

%% GRID GENERATION
xg = makegrid_cart_3D(p);

% these new variables are just for your information, they are written to disk by eq2dist().
[nsi, vs1i, Tsi, xgin, ns, vs1, Ts] = eq2dist(p, xg);

state.nsi = nsi;
state.vs1i = vs1i;

state.Tsi = Tsi;
state.xgin = xgin;
state.ns = ns;
state.vs1 = vs1;
state.Ts = Ts;

%% potential boundary conditions

if p.lxp == 1 || p.lyp == 1
  E = Efield_BCs_2d(p);
else % 3D
  E = Efield_BCs_3d(p);
end

if ~nargout, clear('state', 'E'), end
end % function
