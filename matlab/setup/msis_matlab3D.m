function natm = msis_matlab3D(p, xg)
%% calls MSIS Fortran exectuable from Matlab.
% compiles if not present
%
% [f107a, f107, ap] = activ;
%     COLUMNS OF DATA:
%       1 - ALT
%       2 - HE NUMBER DENSITY(M-3)
%       3 - O NUMBER DENSITY(M-3)
%       4 - N2 NUMBER DENSITY(M-3)
%       5 - O2 NUMBER DENSITY(M-3)
%       6 - AR NUMBER DENSITY(M-3)
%       7 - TOTAL MASS DENSITY(KG/M3)
%       8 - H NUMBER DENSITY(M-3)
%       9 - N NUMBER DENSITY(M-3)
%       10 - Anomalous oxygen NUMBER DENSITY(M-3)
%       11 - TEMPERATURE AT ALT
%
narginchk(2,2)
validateattributes(p, {'struct'}, {'scalar'})
validateattributes(xg, {'struct'}, {'scalar'})

cwd = fileparts(mfilename('fullpath'));
exeloc = [cwd,'/../build'];
makedir(exeloc)
exe = absolute_path([exeloc,'/msis_setup']);
if ispc, exe = [exe, '.exe']; end

if ~is_file(exe)
  src = [absolute_path([cwd,'/../../src/vendor/msis00/msis00_gfortran.f']), ' ', ...
         absolute_path([cwd,'/../../src/neutral/call_msis_gfortran.f90'])];
  % -static avoids problems with missing .so or .dll, from Matlab's
  % internal shell
  fc = getenv('FC');
  if isempty(fc)
    fc = 'gfortran  -static -std=legacy -w';
  end
  cmd = [fc, ' ', src,' -o ',exe];
  disp(cmd)
  status = system(cmd);
  assert(status==0, 'failed to compile MSISe00')
end
assert(is_file(exe), ['MSIS setup executable not found: ', exe])
%% SPECIFY SIZES ETC.
lx1=xg.lx(1); lx2=xg.lx(2); lx3=xg.lx(3);
alt=xg.alt(:)/1e3;
glat=xg.glat(:);
glon=xg.glon(:);
lz=lx1*lx2*lx3;
%% CONVERT DATES/TIMES/INDICES INTO MSIS-FRIENDLY FORMAT
f107a = p.activ(1);
f107 = p.activ(2);
ap = p.activ(3);
ap3 = p.activ(3);
doy = datenum(p.ymd(1), p.ymd(2), p.ymd(3)) - datenum(p.ymd(1),1,1) + 1;

disp(['MSIS00 using DOY:  ',int2str(doy)])
yearshort = mod(p.ymd(1),100);
iyd = yearshort*1000+doy;
%% KLUDGE THE BELOW-ZERO ALTITUDES SO THAT THEY DON'T GIVE INF
alt(alt(:)<=0)=1;
%% FIND A UNIQUE IDENTIFIER FOR THE INPUT FILE
fin = [tempdir, '/msis_setup_input.dat'];
%% CREATE AND INPUT FILE FOR FORTRAN PROGRAM
fid=fopen(fin,'w');
fwrite(fid,iyd,'integer*4');
fwrite(fid,p.UTsec0,'integer*4');
fwrite(fid,f107a,'real*4');
fwrite(fid,f107,'real*4');
fwrite(fid,ap,'real*4');
fwrite(fid,ap3,'real*4');
fwrite(fid,lz,'integer*4');
fwrite(fid,glat,'real*4');
fwrite(fid,glon,'real*4');
fwrite(fid,alt,'real*4');
fclose(fid);
%% CALL MSIS AND READ IN RESULTING BINARY FILE
fout = [tempdir, '/msis_setup_output.dat'];
cmd = [exe,' ',fin,' ',fout,' ',int2str(lz)];
disp(cmd)
[status, msg] = system(cmd);   %output written to file
if status~=0, error(['msis setup failed: ',msg]), end

fid=fopen(fout,'r');
msisdat=fread(fid,lz*11,'real*4=>real*8');
msisdat=reshape(msisdat,[11 lz]);
msisdat=msisdat';
fclose(fid);
%% ORGANIZE
nO=reshape(msisdat(:,3),lx1,lx2,lx3);
nN2=reshape(msisdat(:,4),lx1,lx2,lx3);
nO2=reshape(msisdat(:,5),lx1,lx2,lx3);
Tn=reshape(msisdat(:,11),lx1,lx2,lx3);
nN=reshape(msisdat(:,9),lx1,lx2,lx3);
nNO=4e-1*exp(-3700./Tn).*nO2+5e-7*nO;       %Mitra, 1968
nH=reshape(msisdat(:,8),lx1,lx2,lx3);
natm=cat(4,nO,nN2,nO2,Tn,nN,nNO,nH);

delete(fin);
delete(fout);

if nargout==0, clear('natm'), end
end
