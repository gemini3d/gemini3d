function config()
%% 2D east-west test equilibrium

cwd = fileparts(mfilename('fullpath'));

p.format = 'hdf5';
p.nml = [cwd,'/config.nml'];
p.simdir = [cwd, '/../../../gemini_sim/test2d_eq'];
p.xdist = 1200e3;    %eastward distance
p.ydist = 600e3;    %northward distance
p.lxp = 15;
p.lyp = 1;
p.glat = 67.11;
p.glon = 212.95;
p.I = 90;
p.nmf = 5e11;
p.nme = 2e11;
%% ADD PATHS
gemdir = [cwd, '/../..'];
addpath([gemdir, '/script_utils']);
addpath([gemdir, '/setup']);
addpath([gemdir, '/setup/gridgen'])

p = merge_struct(p, read_config(p.nml));
%% GRID GENERATION
xg = makegrid_cart_3D(p);

%% WRITE GRID & INITIAL CONDITIONS
writegrid(p, xg);

[ns,Ts,vsx1] = eqICs3D(p, xg);
writedata(p.ymd, p.UTsec0, ns, vsx1, Ts, p.simdir, p.format);
end % function
