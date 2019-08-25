%
%function data=msis_matlab(alt,glat,glon,iyd,sec,f107a,f107,ap)
%
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

function natm = msis_matlab3D(xg,UT,dmy,activ)
% [f107a, f107, ap] = activ;
narginchk(4,4)
validateattributes(xg,{'struct'},{'scalar'})
validateattributes(UT,{'numeric'},{'nonnegative','scalar'}, mfilename, "UT decimal hour from midnight", 2)
validateattributes(dmy,{'numeric'},{'positive','vector','numel',3})
validateattributes(activ,{'numeric'},{'positive','vector','numel',3})


cwd = fileparts(mfilename('fullpath'));
exeloc = [cwd,'/../'];
exe = [exeloc,'msis_setup'];

if ~exist(exe,'file'), error('MSIS setup executable not found'), end
%% SPECIFY SIZES ETC.
lx1=xg.lx(1); lx2=xg.lx(2); lx3=xg.lx(3);
alt=xg.alt(:)/1e3;
glat=xg.glat(:);
glon=xg.glon(:);
lz=lx1*lx2*lx3;
%% CONVERT DATES/TIMES/INDICES INTO MSIS-FRIENDLY FORMAT
f107a=activ(1); f107=activ(2); ap=activ(3); ap3=activ(3);
doy = datenum(dmy(3), dmy(2), dmy(1)) - datenum(dmy(3),1,1) + 1;

disp(['MSIS00 using DOY:  ',int2str(doy)])
yearshort = mod(dmy(3),100);
iyd = yearshort*1000+doy;
sec = round(UT*3600);
%% KLUDGE THE BELOW-ZERO ALTITUDES SO THAT THEY DON'T GIVE INF
alt(alt(:)<=0)=1;
%% FIND A UNIQUE IDENTIFIER FOR THE INPUT FILE
fin = tempname;
%% CREATE AND INPUT FILE FOR FORTRAN PROGRAM
fid=fopen(fin,'w');
fwrite(fid,iyd,'integer*4');
fwrite(fid,sec,'integer*4');
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
fout = tempname;
disp(['MSIS00 input: ', fin])
disp(['MSIS00 output: ', fout])

[status, msg] = system([exe,' ',fin,' ',fout]);   %output written to file
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
