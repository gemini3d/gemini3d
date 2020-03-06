function ok = compare_all(outdir, refdir)

% the absolute and relative tolerance account for slight IEEE-754 based differences,
% including non-associativity of floating-point arithmetic.
% these parameters are a bit arbitrary.

% per MZ Oct 17, 2018:
% Ti,Te=1 K
% ne=1e6 m-3
% vi,v2,v3=1 m/s
% J1,J2,J3 = 1e-9

% MZ wants to change what we consider signficant...
% Ti,Te=5 K
% ne=1e7 m-3
% vi,v2,v3=2 m/s
% J1,J2,J3 = 1e-9
narginchk(2,2)

addpath([fileparts(mfilename('fullpath')), '/vis'])

tol.rtol = 1e-5;
tol.rtolN = 1e-5;
tol.rtolT = 1e-5;
tol.rtolJ = 1e-5;
tol.rtolV = 1e-5;
tol.atol=1e-8;
tol.atolN=1e9;
tol.atolT=100;
tol.atolJ=1e-7;
tol.atolV=50;

outdir = absolute_path(outdir);
refdir = absolute_path(refdir);

exist_or_skip(outdir, 'dir')
exist_or_skip(refdir, 'dir')
%% check that paths not the same
if strcmp(outdir, refdir), error([outdir, ' and ', refdir, ' directories are the same']), end
%% READ IN THE SIMULATION INFORMATION
params = read_config([outdir, filesep, 'inputs']);

lxs = simsize(outdir);
lxs_ref = simsize(refdir);
if ~all(lxs==lxs_ref)
  error(['ref dims ',int2str(lxs_ref), ' != this sim dims ', int2str(lxs)])
end

disp(['sim grid dimensions: ',num2str(lxs)])

%% TIMES OF INTEREST
times = params.UTsec0:params.dtout:params.UTsec0 + params.tdur;
Nt = length(times);
assert(Nt > 1, [outdir, ' simulation did not run long enough'])

ymd=params.ymd;
UTsec=params.UTsec0;

ok = false;

for it=1:Nt
  st = ['UTsec ', num2str(times(it))];
  out = loadframe(outdir,ymd,UTsec);
  ref = loadframe(refdir,ymd,UTsec);

  ok = ok + ~assert_allclose(out.ne,ref.ne,tol.rtolN,tol.atolN,['Ne ',st], true);

  if false
    ok = ok + ~assert_allclose(out.v1,ref.v1,tol.rtolV,tol.atolV,['V1 ', st], true);
  end
  ok = ok + ~assert_allclose(out.v2,ref.v2,tol.rtolV,tol.atolV,['V2 ', st], true);
  ok = ok + ~assert_allclose(out.v3,ref.v3,tol.rtolV,tol.atolV,['V3 ', st], true);

  if false
    ok = ok + ~assert_allclose(out.Ti,ref.Ti,tol.rtolT,tol.atolT,['Ti ', st], true);
  end
  ok = ok + ~assert_allclose(out.Te,ref.Te,tol.rtolT,tol.atolT,['Te ', st], true);

  ok = ok + ~assert_allclose(out.J1,ref.J1,tol.rtolJ,tol.atolJ,['J1 ', st], true);
  ok = ok + ~assert_allclose(out.J2,ref.J2,tol.rtolJ,tol.atolJ,['J2 ', st], true);
  ok = ok + ~assert_allclose(out.J3,ref.J3,tol.rtolJ,tol.atolJ,['J3 ', st], true);

  %% assert time steps have unique output (earth always rotating...)
  if it>1
    ok = ok + ~assert_allclose(Ne,out.ne,tol.rtol,tol.atol,['Ne ', st,' too similar to prior step'],true, true);
    %ok = ok + ~assert_allclose(v1,out.v1,tol.rtol,tol.atol,['V1 ', st,' too similar to prior step'],true, true);
    ok = ok + ~assert_allclose(v2,out.v2,tol.rtol,tol.atol,['V2 ', st,' too similar to prior step'],true, true);
    ok = ok + ~assert_allclose(v3,out.v3,tol.rtol,tol.atol,['V3 ', st,' too similar to prior step'],true, true);
  end
  if it==3
   %ok = ok + ~assert_allclose(Ti,out.Ti,tol.rtol,tol.atol,['Ti ', st,' too similar to prior step'],true, true);
    ok = ok + ~assert_allclose(Te,out.Te,tol.rtol,tol.atol,['Te ', st,' too similar to prior step'],true, true);
  end
  if it==2
    ok = ok + ~assert_allclose(J1,out.J1,tol.rtol,tol.atol,['J1 ', st,' too similar to prior step'],true,true, true);
    ok = ok + ~assert_allclose(J2,out.J2,tol.rtol,tol.atol,['J2 ', st,' too similar to prior step'],true,true, true);
    ok = ok + ~assert_allclose(J3,out.J3,tol.rtol,tol.atol,['J3 ', st,' too similar to prior step'],true,true, true);
  end

  Ne = out.ne; v1=out.v1; v2=out.v2; v3=out.v3; Ti=out.Ti; Te=out.Te; J1=out.J1; J2=out.J2; J3=out.J3;

  [ymd,UTsec] = dateinc(params.dtout,ymd,UTsec);

end

if ok ~= 0
  error([int2str(ok), ' compare errors'])
else
  disp(['OK: Gemini output comparison of ',int2str(Nt),' time steps.'])
end

if nargout==0, clear('ok'), end

end % function
