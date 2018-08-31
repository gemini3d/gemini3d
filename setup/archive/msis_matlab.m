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

function natm=msis_matlab(xg,UT,dmy,activ)

exeloc='./neutral/';
lx1=xg.lx(1); lx2=xg.lx(2);
alt=xg.alt(:)/1e3;
glat=xg.glat(:);
glon=xg.glon(:);
lz=lx1*lx2;


%CONVERT DATES/TIMES/INDICES INTO MSIS-FRIENDLY FORMAT
f107a=activ(1); f107=activ(2); ap=activ(3); ap3=activ(3);
dom=dmy(1); month=dmy(2); year=dmy(3);
%doy=round(30.3635*(month-1)+dom);
monthdays=[31,28,31,30,31,30,31,31,30,31,30,31];
if mod(year,4)==0
  monthdays(2)=29;    %leap year!
end
if month==1
  doy=dom;
else
  doy=sum(monthdays(1:month-1))+dom;
end
fprintf('\nMSIS00 using DOY:  %d \n',doy);
yearshort=mod(year,100);
iyd=yearshort*1000+doy;
t0=UT*3600;
sec=round(t0);


%KLUDGE THE BELOW-ZERO ALTITUDES SO THAT THEY DON'T GIVE INF
inds=find(alt(:)<=0);
alt(inds)=1;


%FIND A UNIQUE IDENTIFIER FOR THE INPUT FILE
if exist('getpid')
   pid=getpid;
else
   pid=feature('getpid');
end
filename=['msisinput_',num2str(pid),'.dat'];


%CREATE AND INPUT FILE FOR FORTRAN PROGRAM (BINARY WOULD BE FASTER!!!)
fid=fopen([exeloc,filename],'w');
fprintf(fid,'%d %d\n',iyd,sec);
fprintf(fid,'%f %f\n',f107a,f107);
fprintf(fid,'%f %f\n',ap,ap3);
fprintf(fid,'%d\n',lz);
for iz=1:lz
    fprintf(fid,'%f %f %f\n',glat(iz),glon(iz),alt(iz));
end
fclose(fid);


%MAKE THE CALL TO THE SHELL FOR THE DATA FOR THESE ALTS
if ismac
  path1 = getenv('PATH');
  path1 = [path1 ':/usr/local/bin'];
  setenv('PATH', path1);
  setenv('DYLD_LIBRARY_PATH', '/usr/local/bin')
end
%[flg,datchar]=system([exeloc,'msis ',exeloc,filename]);    %output written to shell env.
%msisdat=str2num(datchar);


%CALL MSIS AND READ IN RESULTING BINARY FILE
fprintf('MSIS00 using binary file output');
filename2=['msisoutput_',num2str(pid),'.dat'];
[flg,datchar]=system([exeloc,'msis ',exeloc,filename,' ',exeloc,filename2]);   %output written to file
fid=fopen([exeloc,filename2],'r');
msisdat=fread(fid,lz*11,'real*4=>real*8');
lz,11*lz
size(msisdat)
msisdat=reshape(msisdat,[11 lz]);
msisdat=msisdat';
fclose(fid);

%GARBAGE CLEANUP
system(['rm ',exeloc,filename]);
system(['rm ',exeloc,filename2]);


%ORGANIZE
nO=reshape(msisdat(:,3),lx1,lx2);
nN2=reshape(msisdat(:,4),lx1,lx2);
nO2=reshape(msisdat(:,5),lx1,lx2);
Tn=reshape(msisdat(:,11),lx1,lx2);
nN=reshape(msisdat(:,9),lx1,lx2);
nNO=4e-1*exp(-3700./Tn).*nO2+5e-7*nO;       %Mitra, 1968
nH=reshape(msisdat(:,8),lx1,lx2);
natm=cat(3,nO,nN2,nO2,Tn,nN,nNO,nH);

end
