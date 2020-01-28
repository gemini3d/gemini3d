%% LOWRES east-west 2D EXAMPLE FOR TESTING

cwd = fileparts(mfilename('fullpath'));

p.format = 'hdf5';
p.simID = 'test2d';
p.nml = [cwd,'/config.nml'];
p.xdist = 200e3;    %eastward distance
p.ydist = 600e3;    %northward distance
p.lxp = 80;
p.lyp = 1;
p.glat = 67.11;
p.glon = 212.95;
p.I = 90;
eq_dir = '../../../gemini_sim/test2d_fang';

%% ADD PATHS
gemdir = [cwd, '/../..'];
addpath([gemdir, '/setup']);
addpath([gemdir, '/setup/gridgen'])

cfg = read_config(p.nml);
UThour = cfg.UTsec0 / 3600;
dmy = flip(cfg.ymd(:).');
% import geomagindices as gi
% gi.get_indices('2013-02-20T05', smoothdays=81)
activ = [108.9, 111.0, 5];
nmf=5e11;
nme=2e11;
%% GRID GENERATION
xg = makegrid_cart_3D(p.xdist, p.lxp, p.ydist, p.lyp, p.I, p.glat, p.glon);
writegrid(xg, [eq_dir, '/inputs'], p.format);
%% GENERATE INITIAL CONDITIONS FOR A PARTICULAR EVENT
[ns,Ts,vsx1] = eqICs3D(xg,UThour,dmy,activ,nmf,nme);

writedata(dmy, cfg.UTsec0, ns, vsx1, Ts, [eq_dir, '/inputs'], p.format)
