function [nsi,vs1i,Tsi,xgin,ns,vs1,Ts] = eq2dist(eqdir, simID, xg, format, outdir)

narginchk(4, 5)
validateattributes(eqdir, {'char', 'string'}, {'vector'})
validateattributes(simID, {'char', 'string'}, {'vector'})
validateattributes(xg, {'struct'}, {'scalar'}, mfilename, 'grid struct',3)
validateattributes(format, {'char'}, {'vector'}, mfilename, 'raw or hdf5',4)
if nargin < 5, outdir = []; end
%% READ SIMULATION INFORMATION
[ymd0,UTsec0,tdur,dtout,flagoutput,mloc] = readconfig([eqdir, '/inputs']);
xgin = readgrid([eqdir, '/inputs']);

%% FIND THE DATE OF THE END FRAME OF THE SIMULATION
% PRESUMABLY THIS WILL BE THE STARTING point FOR another
[ymdend,UTsecend] = dateinc(tdur,ymd0,UTsec0);

%% LOAD THE FRAME
[ne,mlatsrc,mlonsrc,xgin,v1,Ti,Te,J1,v2,v3,J2,J3,filename,Phitop,ns,vs1,Ts] = ...
    loadframe(eqdir,ymdend,UTsecend,flagoutput,mloc,xgin);

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
% this uses SIMID as output directory and filename tag
if isempty(outdir)
  basedir = [eqdir,'/../input/'];
  outdir = [basedir, simID];
end
writegrid(xg, outdir, format);
dmy=[ymdend(3),ymdend(2),ymdend(1)];
writedata(dmy,UTsecend,nsi,vs1i,Tsi,outdir,simID, format);

end % function eq2dist
