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

function natm=msis_matlab3D(xg,UT,dmy,activ)

%exeloc='./';
dat=which('msis_matlab3D');
rest=dat;
piece=[];
strpath=[];
while ~isempty(rest)
    strpath=[strpath,'/',piece];
    [piece,rest]=strtok(rest,'/');
end
exeloc=[strpath,'/'];


%SPECIFY SIZES ETC.
lx1=xg.lx(1); lx2=xg.lx(2); lx3=xg.lx(3);
alt=xg.alt(:)/1e3;
glat=xg.glat(:);
glon=xg.glon(:);
lz=lx1*lx2*lx3;


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


%CREATE AND INPUT FILE FOR FORTRAN PROGRAM
fid=fopen([exeloc,filename],'w');
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

%{
fprintf(fid,'%d %d\n',iyd,sec);
fprintf(fid,'%f %f\n',f107a,f107);
fprintf(fid,'%f %f\n',ap,ap3);
fprintf(fid,'%d\n',lz);
for iz=1:lz
    fprintf(fid,'%f %f %f\n',glat(iz),glon(iz),alt(iz));
end
%}
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
system_command=[exeloc,'msis ',exeloc,filename,' ',exeloc,filename2];
[flg,datchar]=system(system_command);   %output written to file
fid=fopen([exeloc,filename2],'r');
msisdat=fread(fid,lz*11,'real*4=>real*8');
msisdat=reshape(msisdat,[11 lz]);
msisdat=msisdat';
fclose(fid);


%GARBAGE CLEANUP
system(['rm ',exeloc,filename]);
system(['rm ',exeloc,filename2]);


%ORGANIZE
nO=reshape(msisdat(:,3),lx1,lx2,lx3);
nN2=reshape(msisdat(:,4),lx1,lx2,lx3);
nO2=reshape(msisdat(:,5),lx1,lx2,lx3);
Tn=reshape(msisdat(:,11),lx1,lx2,lx3);
nN=reshape(msisdat(:,9),lx1,lx2,lx3);
nNO=4e-1*exp(-3700./Tn).*nO2+5e-7*nO;       %Mitra, 1968
nH=reshape(msisdat(:,8),lx1,lx2,lx3);
natm=cat(4,nO,nN2,nO2,Tn,nN,nNO,nH);

end
