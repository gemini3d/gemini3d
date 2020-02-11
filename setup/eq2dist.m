function [nsi,vs1i,Tsi,xgin,ns,vs1,Ts] = eq2dist(p, xg)
% read and interpolate equilibrium simulation data, writing new
% interpolated grid.
narginchk(2, 2)
validateattributes(p, {'struct'}, {'scalar'}, mfilename, 'parameters', 1)
validateattributes(xg, {'struct'}, {'scalar'}, mfilename, 'grid struct', 2)
%% Paths
% this script is called from numerous places, so ensure necessary path
cwd = fileparts(mfilename('fullpath'));
addpath([cwd,'/../vis'])
%% READ Equilibrium SIMULATION INFO
[ymd0,UTsec0,tdur,dtout, flagoutput, mloc] = readconfig(p.eqnml);
xgin = readgrid([p.eqdir, '/inputs'], p.format);

%% FIND THE DATE OF THE END FRAME OF THE SIMULATION
% PRESUMABLY THIS WILL BE THE STARTING point FOR another
[ymdend,UTsecend] = dateinc(tdur,ymd0,UTsec0);

%% LOAD THE FRAME
[ne,mlatsrc,mlonsrc,xgin,v1,Ti,Te,J1,v2,v3,J2,J3,filename,Phitop,ns,vs1,Ts] = ...
    loadframe(p.eqdir, ymdend, UTsecend, flagoutput, mloc, xgin, [], p.eqnml);

%% check input to interpolation
assert(all(isfinite(ns(:))), 'non-finite density')
assert(all(isfinite(vs1(:))), 'non-finite drift')
assert(all(isfinite(Ts(:))), 'non-finite temperature')

%% DO THE INTERPOLATION
[nsi,vs1i,Tsi] = model_resample(xgin, ns, vs1, Ts, xg);

%% check IF THE INTERPOLATION WENT WEIRD...
assert(all(isfinite(nsi(:))), 'non-finite interpolated density')
assert(all(isfinite(vs1i(:))), 'non-finite interpolated drift')
assert(all(isfinite(Tsi(:))), 'non-finite interpolated temperature')

%% WRITE OUT THE GRID
writegrid(xg, p.simdir, p.format);
dmy=[ymdend(3),ymdend(2),ymdend(1)];
writedata(dmy,UTsecend,nsi,vs1i,Tsi, p.simdir, p.format);

end % function eq2dist
