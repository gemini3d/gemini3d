function config()
%% 2D east-west test
%
cwd = fileparts(mfilename('fullpath'));

p.format = 'hdf5';
p.nml = [cwd,'/config.nml'];
p.eqdir = [cwd, '/../../../gemini_sim/test2d_eq'];
p.xdist = 200e3;    %eastward distance
p.ydist = 600e3;    %northward distance
p.lxp = 80;
p.lyp = 1;
p.glat = 67.11;
p.glon = 212.95;
p.I = 90;

%% ADD PATHS
gemdir = [cwd, '/../..'];
addpath([gemdir, '/script_utils']);
addpath([gemdir, '/setup']);
addpath([gemdir, '/setup/gridgen'])

p = merge_struct(p, read_config(p.nml));
% import geomagindices as gi
% gi.get_indices('2013-02-20T05', smoothdays=81)
% p.activ = [108.9, 111.0, 5];

%% GRID GENERATION
xg = makegrid_cart_3D(p);

%% WRITE GRID & INITIAL CONDITIONS
writegrid(p, xg);

[ns,Ts,vsx1] = eqICs3D(p, xg);
writedata(p.ymd, p.UTsec0, ns, vsx1, Ts, p.simdir, p.format);
end % function
