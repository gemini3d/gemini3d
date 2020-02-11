%% LOWRES 3D EXAMPLE FOR TESTING

cwd = fileparts(mfilename('fullpath'));

p.format = 'hdf5';
p.nml = [cwd,'/config.nml'];
p.eqdir = [cwd, '../../gemini_sim/test3d_eq'];
p.xdist = 200e3;    %eastward distance
p.ydist = 200e3;    %northward distance
vlxp = 20;
p.lyp = 20;
p.glat = 67.11;
p.glon = 212.95;
p.I = 90;

%% ADD PATHS
gemdir = [cwd, '/../..'];
addpath([gemdir, '/setup']);
addpath([gemdir, '/setup/gridgen'])


%% GRID GENERATION
xg = makegrid_cart_3D_lowresx1(p);

%% IC FROM EQ SIMULATION
[nsi,vs1i,Tsi,xgin,ns,vs1,Ts] = eq2dist(p, xg);
