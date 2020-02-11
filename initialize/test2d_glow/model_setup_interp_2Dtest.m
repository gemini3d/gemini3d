%% LOWRES east-west 2D EXAMPLE FOR TESTING

cwd = fileparts(mfilename('fullpath'));

p.format = 'hdf5';
p.nml = [cwd,'/config.nml'];
p.eqdir = [cwd, '../../gemini_sim/test3d_eq'];
p.xdist = 200e3;    %eastward distance
p.ydist = 600e3;    %northward distance
p.lxp = 80;
p.lyp = 1;
p.glat = 67.11;
p.glon = 212.95;
p.I = 90;

%% ADD PATHS
gemdir = [cwd, '/../..'];
addpath([gemdir, '/setup']);
addpath([gemdir, '/setup/gridgen'])


%% GRID GENERATION
xg = makegrid_cart_3D(p);

%% IC FROM EQ SIMULATION
[nsi,vs1i,Tsi,xgin,ns,vs1,Ts] = eq2dist(p, xg);
