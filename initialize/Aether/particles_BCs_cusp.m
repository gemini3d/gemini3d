%addpath ./restore_idl;
addpath ../../script_utils;


%REFERENCE GRID TO USE
direcconfig='./'
direcgrid='~/zettergmdata/simulations/input/Aether/'

%CREATE SOME SPACE FOR OUTPUT FILES
outdir='~/zettergmdata/simulations/input/Aether_particles_cusp/';
system(['mkdir ',outdir]);
system(['rm ',outdir,'/*']);


%READ IN THE SIMULATION INFORMATION (MEANS WE NEED TO CREATE THIS FOR THE SIMULATION WE WANT TO DO)
if (~exist('ymd0','var'))
  [ymd0,UTsec0,tdur,dtout,flagoutput,mloc]=readconfig([direcconfig,'/config.ini']);
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
mlonsig=1/4*(mlonmax-mlonmin);



%TIME VARIABLE (SECONDS FROM SIMULATION BEGINNING)
tmin=0;
tmax=tdur;
%lt=tdur+1;
%time=linspace(tmin,tmax,lt)';
time=tmin:1:tmax;
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


%CREATE THE PRECIPITAITON INPUT DATA
Qit=zeros(llon,llat,lt);
E0it=zeros(llon,llat,lt);

mlonsig=1.5;
mlatsig=0.25;

Qpk=2;
E0pk=100;
for it=1:lt
  shapefn=exp(-(MLON-mlonmean).^8/2/mlonsig^8).*exp(-(MLAT-mlatmean).^2/2/mlatsig^2);
  Qittmp=Qpk.*shapefn;
  E0it(:,:,it)=E0pk;%*ones(llon,llat);     %eV
  inds=find(Qittmp<0.05);    %define a background flux (enforces a floor for production rates)
  Qittmp(inds)=0.05;
  Qit(:,:,it)=Qittmp;
end


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


%RESTORE PATH
%rmpath ./restore_idl;
rmpath ../../script_utils;
