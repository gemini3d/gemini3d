function model_setup_equilibrium(p)
narginchk(1,1)
validateattributes(p, {'struct'}, {'scalar'}, mfilename, 'parameters', 1)

%% ADD PATHS FOR FUNCTIONS
cwd = fileparts(mfilename('fullpath'));
gemdir = [cwd, '/../../gemini'];
assert(isfolder(gemdir), [gemdir, ' not found'])
for d = {'script_utils', 'setup', 'setup/gridgen', 'vis'}
  addpath([gemdir, filesep, d{:}])
end

%% GRID GENERATION
xg = makegrid_cart_3D(p);
writegrid(p, xg);

%% Equilibrium input generation
cfg = read_config(p.nml);
[ns,Ts,vsx1] = eqICs3D(p, xg);
% Note: should be rewritten to include the neutral module form the fortran code
writedata(p.ymd, cfg.UTsec0, ns, vsx1, Ts, p.simdir, p.format);

end % function
