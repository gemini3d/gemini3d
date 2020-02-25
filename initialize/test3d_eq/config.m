function config()
%% LOWRES 3D EXAMPLE FOR TESTING
%

%% boilerplate for each config.m file
cwd = fileparts(mfilename('fullpath'));
gemroot = getenv('GEMINI_ROOT');
addpath([gemroot, '/script_utils'], [gemroot,'/setup'])
p.nml = [cwd,'/config.nml'];
p = merge_struct(p, read_nml(p.nml));
p.simdir = absolute_path(['../', fileparts(p.indat_size)]);

p.realbits = 64;
p.format = 'h5';

%% setup simulation
model_setup_equilibrium(p)

end % function
