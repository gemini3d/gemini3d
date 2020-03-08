function [nsi,vs1i,Tsi] = eq2dist(p, xg)
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
peq = read_config(p.eqdir);
xgin = readgrid(p.eqdir, p.format, p.realbits);


%% END FRAME time of equilibrium simulation
% PRESUMABLY THIS WILL BE THE STARTING point FOR another
[ymd_end,UTsec_end] = dateinc(peq.tdur,peq.ymd,peq.UTsec0);


%% LOAD THE last equilibrium frame
dat = loadframe(p.eqdir, ymd_end, UTsec_end, peq.flagoutput, peq.mloc, xgin, p.format, p.eqdir);


%% sanity check equilibrium simulation input to interpolation
check_density(dat.ns)
check_drift(dat.vs1)
check_temperature(dat.Ts)


%% DO THE INTERPOLATION
[nsi,vs1i,Tsi] = model_resample(xgin, dat.ns, dat.vs1, dat.Ts, xg);


%% sanity check interpolated variables
check_density(nsi)
check_drift(vs1i)
check_temperature(Tsi)


%% WRITE OUT THE GRID
writegrid(p, xg)
writedata(ymd_end, UTsec_end,nsi,vs1i,Tsi, p.simdir, p.format, p.realbits);

end % function eq2dist


function check_density(n)
narginchk(1,1)

n = n(:);
assert(all(isfinite(n)), 'non-finite density')
assert(all(n > 0), 'negative density')
assert(max(n > 1e6), 'too small maximum density')
end

function check_drift(v)
narginchk(1,1)

v = v(:);
assert(all(isfinite(v)), 'non-finite drift')
assert(all(abs(v) < 10e3), 'excessive drift velocity')
end

function check_temperature(T)

T = T(:);
assert(all(isfinite(T)), 'non-finite temperature')
assert(all(T > 0), 'negative temperature')
assert(max(T) > 500, 'too cold maximum temperature')

end
