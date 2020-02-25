function config()
%% 3D test
%

%% boilerplate for each config.m file
cwd = fileparts(mfilename('fullpath'));
gemroot = getenv('GEMINI_ROOT');
addpath([gemroot, '/script_utils'], [gemroot,'/setup'])
p.nml = [cwd,'/config.nml'];
p = merge_struct(p, read_nml(p.nml));
p.simdir = absolute_path(['../', fileparts(p.indat_size)]);

p.eqdir = absolute_path([p.simdir, '/../test3d_eq']);
p.eqnml = [p.eqdir, '/inputs/config.nml'];

p.format = 'h5';
p.realbits = 64;
%% perturbations
p.Etarg = 50e-3;  % [V/m]

% normalized widths
p.Efield_fracwidth = 1/7;
p.precip_latwidth = 1/4;
p.precip_lonwidth = 1/4;

%% setup simulation
model_setup_interp(p)

end % function
