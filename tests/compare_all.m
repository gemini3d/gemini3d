function compare_all(dir1, dir2)
  
  % the absolute and relative tolerance account for slight IEEE-754 based differences,
  % including non-associativity of floating-point arithmetic.
  % these parameters are a bit arbitrary.
  
  % ==== Haswell Ubuntu 18.04 vs. Matt's March 20,2018 test2d
%   rtol=1e-8; rtolN=0.0001; rtolT=0.0001; rtolJ=0.01;   rtolV=0.1; 
%   atol=1e-8;                             atolJ=0.0001; atolV=0.0001;    
 
  % ==== phys_const.f90 update f225127
%   rtol=1e-7; rtolN=0.001; rtolT=0.001; rtolJ=rtol;   rtolV=0.05; 
%   atol=1e-7;                             atolJ=atol; atolV=0.1;  
 
  % ==== Haswell Ubuntu 18.04 vs. itself on various git revs
  
narginchk(2,2)
validateattributes(dir1, {'char'}, {'vector'}, mfilename,'directory to compare',1)
validateattributes(dir2, {'char'}, {'vector'}, mfilename,'directory to compare',2)
  
rtol=1e-6; rtolN=1e-3; rtolT=1e-4; rtolJ=0.1; rtolV=rtol; 
atol=1e-9;                        atolJ=atol; atolV=atol;    

cwd = fileparts(mfilename('fullpath'));
    
addpath([cwd, filesep, '..', filesep, 'script_utils'])
addpath([cwd, filesep, '..', filesep, 'vis'])
%% READ IN THE SIMULATION INFORMATION
[ymd0,UTsec0,tdur,dtout,~,mloc]=readconfig([dir1,filesep,'inputs/config.ini']);
%% load grid
xg=readgrid([dir1,filesep,'inputs',filesep]);
%% TIMES OF INTEREST
times=UTsec0:dtout:UTsec0+tdur;
Nt = length(times);
assert(Nt > 1, 'simulation did not run long enough')

ymd=ymd0;
UTsec=UTsec0;
  
for it=1:Nt 
  [neA,v1A,TiA,TeA,J1A,v2A,v3A,J2A,J3A] = loadframe(dir1,UTsec,ymd,UTsec0,ymd0,mloc,xg);
  [neB,v1B,TiB,TeB,J1B,v2B,v3B,J2B,J3B] = loadframe(dir2,UTsec,ymd,UTsec0,ymd0,mloc,xg);
  
  assert_allclose(neA,neB,rtolN,atol,['Ne accuracy, time #',int2str(it)])
  
  if false
    assert_allclose(v1A,v1B,rtolV,atolV,['V1 accuracy, time #',int2str(it)])
  end
  assert_allclose(v2A,v2B,rtolV,atolV,['V2 accuracy, time #',int2str(it)])
  assert_allclose(v3A,v3B,rtolV,atolV,['V3 accuracy, time #',int2str(it)])
  
  if false
    assert_allclose(TiA,TiB,rtolT,atol,['Ti accuracy, time #',int2str(it)])
  end
  assert_allclose(TeA,TeB,rtolT,atol,['Te accuracy, time #',int2str(it)])
  
  assert_allclose(J1A,J1B,rtolJ,atolJ,['J1 accuracy, time #',int2str(it)])
  assert_allclose(J2A,J2B,rtolJ,atolJ,['J2 accuracy, time #',int2str(it)])
  assert_allclose(J3A,J3B,rtolJ,atolJ,['J3 accuracy, time #',int2str(it)])
  
  %% assert time steps have unique output (earth always rotating...)
  if it>1
    assert_allclose(Ne,neA,rtol,atol,['Ne time step ',int2str(it),' too similar to prior step'],true) 
    %assert_allclose(v1,v1A,rtol,atol,['V1 time step ',int2str(it),' too similar to prior step'],true) 
    assert_allclose(v2,v2A,rtol,atol,['V2 time step ',int2str(it),' too similar to prior step'],true) 
    assert_allclose(v3,v3A,rtol,atol,['V3 time step ',int2str(it),' too similar to prior step'],true) 
  end
  if it==3
   % assert_allclose(Ti,TiA,rtol,atol,['Ti time step ',int2str(it),' too similar to prior step'],true) 
    assert_allclose(Te,TeA,rtol,atol,['Te time step ',int2str(it),' too similar to prior step'],true) 
  end
  if it==2
    assert_allclose(J1,J1A,rtol,atol,['J1 time step ',int2str(it),' too similar to prior step'],true,true) 
    assert_allclose(J2,J2A,rtol,atol,['J2 time step ',int2str(it),' too similar to prior step'],true,true) 
    assert_allclose(J3,J3A,rtol,atol,['J3 time step ',int2str(it),' too similar to prior step'],true,true) 
  end
  
  Ne = neA; v1=v1A; v2=v2A; v3=v3A; Ti=TiA; Te=TeA; J1=J1A; J2=J2A; J3=J3A;
  
  [ymd,UTsec]=dateinc(dtout,ymd,UTsec);
end

disp(['OK: Gemini output comparison of ',int2str(Nt),' time steps.'])

end % function

