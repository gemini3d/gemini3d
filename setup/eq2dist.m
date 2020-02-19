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
xgin = readgrid(p.eqnml, p.format);

%% FIND THE DATE OF THE END FRAME OF THE SIMULATION
% PRESUMABLY THIS WILL BE THE STARTING point FOR another
[ymd_end,UTsec_end] = dateinc(tdur,ymd0,UTsec0);

%% LOAD THE FRAME
dat = loadframe(p.eqdir, ymd_end, UTsec_end, flagoutput, mloc, xgin, [], p.eqnml);

%% check input to interpolation
assert(all(isfinite(dat.ns(:))), 'non-finite density')
assert(all(isfinite(dat.vs1(:))), 'non-finite drift')
assert(all(isfinite(dat.Ts(:))), 'non-finite temperature')

%% DO THE INTERPOLATION
[nsi,vs1i,Tsi] = model_resample(xgin, dat.ns, dat.vs1, dat.Ts, xg);

%% check IF THE INTERPOLATION WENT WEIRD...
assert(all(isfinite(nsi(:))), 'non-finite interpolated density')
assert(all(isfinite(vs1i(:))), 'non-finite interpolated drift')
assert(all(isfinite(Tsi(:))), 'non-finite interpolated temperature')

%% WRITE OUT THE GRID
writegrid(p, xg)
writedata(ymd_end, UTsec_end,nsi,vs1i,Tsi, p.simdir, p.format);

end % function eq2dist
