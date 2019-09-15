%addpath ./restore_idl;
addpath ../../script_utils;


%REFERENCE GRID TO USE
direcconfig='./'
direcgrid='../tests/data/zenodo2d/'

%CREATE SOME SPACE FOR OUTPUT FILES
outdir='../simulations/input/test2d_particles/';
system(['mkdir ',outdir]);
system(['rm ',outdir,'/*']);


%{
%READ IN THE IDL SAVE FILE
%fname='isinglass_eflux_asi.sav';
%fname='isinglass_eflux_asi_sync.sav';
fname='isinglass_eflux_asi_resync.sav';
outargs=restore_idl(fname);
time=outargs.NEW_TIME;
lat=outargs.NEW_LAT;
lon=outargs.NEW_LON;
Qdat=outargs.RESAMP_Q;
E0dat=outargs.RESAMP_EO;
[lt,llon,llat]=size(Qdat);
%}


%READ IN THE SIMULATION INFORMATION (MEANS WE NEED TO CREATE THIS FOR THE SIMULATION WE WANT TO DO)
if (~exist('ymd0','var'))
  [ymd0,UTsec0,tdur,dtout,flagoutput,mloc]=readconfig(direcconfig);
  fprintf('Input config.dat file loaded.\n');
end


%CHECK WHETHER WE NEED TO RELOAD THE GRID (SO THIS ALREADY NEEDS TO BE MADE, AS WELL)
if (~exist('xg','var'))
  %WE ALSO NEED TO LOAD THE GRID FILE
  xg=readgrid([direcgrid,'/']);
  fprintf('Grid loaded.\n');
end


%CREATE A 'DATASET' OF PRECIPITATION CHARACTERISTICS
llon=100;
llat=100;
if (xg.lx(2)==1)    %this is cartesian-specific code
    llon=1;
elseif (xg.lx(3)==1)
    llat=1;
end
thetamin=min(xg.theta(:));
thetamax=max(xg.theta(:));
mlatmin=90-thetamax*180/pi;
mlatmax=90-thetamin*180/pi;
mlonmin=min(xg.phi(:))*180/pi;
mlonmax=max(xg.phi(:))*180/pi;
%mlat=linspace(mlatmin,mlatmax,llat);
%mlon=linspace(mlonmin,mlonmax,llon);
latbuf=1/100*(mlatmax-mlatmin);
lonbuf=1/100*(mlonmax-mlonmin);
mlat=linspace(mlatmin-latbuf,mlatmax+latbuf,llat);
mlon=linspace(mlonmin-lonbuf,mlonmax+lonbuf,llon);
[MLON,MLAT]=ndgrid(mlon,mlat);
mlonmean=mean(mlon);
mlatmean=mean(mlat);



%WIDTH OF THE DISTURBANCE
mlatsig=1/4*(mlatmax-mlatmin);
mlatsig=max(mlatsig,0.01);    %can't let this go to zero...
mlonsig=1/4*(mlonmax-mlonmin);



%TIME VARIABLE (SECONDS FROM SIMULATION BEGINNING)
tmin=0;
tmax=tdur;
%lt=tdur+1;
%time=linspace(tmin,tmax,lt)';
time=tmin:5:tmax;
lt=numel(time);


%COMPUTE THE TOTAL ENERGY FLUX AND CHAR. ENERGY
%{
Q=zeros(lt,llon,llat);
E0=zeros(lt,llon,llat);
for it=1:lt
  Qtmp=squeeze(Qdat(it,:,:));
  Qtmp=interp2(glat,gloncorrected,Qtmp,GLAT(:),GLON(:));    %data are plaid in geographic so do the interpolation in that variable
  Q(it,:,:)=reshape(Qtmp,[1 llon llat]);
  E0tmp=squeeze(E0dat(it,:,:));
  E0tmp=interp2(glat,gloncorrected,E0tmp,GLAT(:),GLON(:));    %data are plaid in geographic so do the interpolation in that variable
  E0(it,:,:)=reshape(E0tmp,[1 llon llat]);
end
%}


%SET UP TIME VARIABLES
ymd=ymd0;
UTsec=UTsec0+time;     %time given in file is the seconds from beginning of hour
UThrs=UTsec/3600;
expdate=cat(2,repmat(ymd,[lt,1]),UThrs(:),zeros(lt,1),zeros(lt,1));
t=datenum(expdate);

%{
%INTERPOLATE/DECIMATE TO 1 SECOND RESOLUTION
samplerate=1;    %sampling rate in seconds
startdate=[ymd0,UTsec0/3600,0,0];
startt=datenum(startdate);
%tmin=ceil(min(t)*86400)/86400;
tmin=ceil(startt*86400)/86400;
tmax=floor(max(t)*86400)/86400;
outputt=tmin:samplerate/86400:tmax;
outputdate=datevec(outputt);
ltout=numel(outputt);
Qit=zeros(llon,llat,ltout);
E0it=zeros(llon,llat,ltout);
for ilat=1:llat
   for ilon=1:llon
       Qhere=squeeze(Qsmooth(:,ilon,ilat));
       E0here=squeeze(E0smooth(:,ilon,ilat));
       Qi=interp1(t,Qhere,outputt);
       E0i=interp1(t,E0here,outputt);
       Qit(ilon,ilat,:)=reshape(Qi,[1 1 ltout]);
       E0it(ilon,ilat,:)=reshape(E0i,[1 1 ltout]);
   end
end
%}


%CREATE THE PRECIPITAITON INPUT DATA
Qit=zeros(llon,llat,lt);
E0it=zeros(llon,llat,lt);
for it=1:lt
   Qit(:,:,it)=10*exp(-(MLON-mlonmean).^2/(2*mlonsig^2)).*exp(-(MLAT-mlatmean).^2/(2*mlatsig^2));         %mW/m^2
%  Qit(:,:,it)=5;
  E0it(:,:,it)=5e3;%*ones(llon,llat);     %eV
end



%%CONVER THE ENERGY TO EV
%E0it=max(E0it,0.100);
%E0it=E0it*1e3;



%SAVE THIS DATA TO APPROPRIATE FILES - LEAVE THE SPATIAL AND TEMPORAL INTERPOLATION TO THE
%FORTRAN CODE IN CASE DIFFERENT GRIDS NEED TO BE TRIED.  THE EFIELD DATA DO
%NO NEED TO BE SMOOTHED.
filename=[outdir,'simsize.dat'];
fid=fopen(filename,'w');
fwrite(fid,llon,'integer*4');
fwrite(fid,llat,'integer*4');
fclose(fid);
filename=[outdir,'simgrid.dat'];
fid=fopen(filename,'w');
fwrite(fid,mlon,'real*8');
fwrite(fid,mlat,'real*8');
fclose(fid);
for it=1:lt
    UTsec=expdate(it,4)*3600+expdate(it,5)*60+expdate(it,6);
    ymd=expdate(it,1:3);
    filename=datelab(ymd,UTsec);
    filename=[outdir,filename,'.dat']
    fid=fopen(filename,'w');
    fwrite(fid,Qit(:,:,it),'real*8');
    fwrite(fid,E0it(:,:,it),'real*8');
    fclose(fid);
end


%ALSO SAVE TO A  MATLAB FILE
save([outdir,'particles.mat'],'mlon','mlat','Qit','E0it','expdate');
