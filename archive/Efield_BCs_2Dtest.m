%REFERENCE GRID TO USE
direcconfig='./'
direcgrid='../tests/data/zenodo2d/input/'


%OUTPUT FILE LOCATION
outdir='../simulations/input/test2d_fields/';
mkdir(outdir);


%READ IN THE SIMULATION INFORMATION (MEANS WE NEED TO CREATE THIS FOR THE SIMULATION WE WANT TO DO)
if (~exist('ymd0','var'))
  [ymd0,UTsec0,tdur,dtout,flagoutput,mloc]=readconfig(direcconfig);
  fprintf('Input config.dat file loaded.\n');
end


%CHECK WHETHER WE NEED TO RELOAD THE GRID (SO THIS ALREADY NEEDS TO BE MADE, AS WELL)
if (~exist('xg','var'))
  %WE ALSO NEED TO LOAD THE GRID FILE
  xg=readgrid([direcgrid,'/']);
  lx1=xg.lx(1); lx2=xg.lx(2); lx3=xg.lx(3);
  fprintf('Grid loaded.\n');
end


%CREATE A 'DATASET' OF ELECTRIC FIELD INFO
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
latbuf=1/100*(mlatmax-mlatmin);
lonbuf=1/100*(mlonmax-mlonmin);
mlat=linspace(mlatmin-latbuf,mlatmax+latbuf,llat);
mlon=linspace(mlonmin-lonbuf,mlonmax+lonbuf,llon);
[MLON,MLAT]=ndgrid(mlon,mlat);
mlonmean=mean(mlon);
mlatmean=mean(mlat);



%WIDTH OF THE DISTURBANCE
fracwidth=1/7;
mlatsig=fracwidth*(mlatmax-mlatmin);
mlonsig=fracwidth*(mlonmax-mlonmin);
sigx2=fracwidth*(max(xg.x2)-min(xg.x2));


%TIME VARIABLE (SECONDS FROM SIMULATION BEGINNING)
tmin=0;
tmax=300;
lt=301;
time=linspace(tmin,tmax,lt)';


%SET UP TIME VARIABLES
ymd=ymd0;
UTsec=UTsec0+time;     %time given in file is the seconds from beginning of hour
UThrs=UTsec/3600;
expdate=cat(2,repmat(ymd,[lt,1]),UThrs,zeros(lt,1),zeros(lt,1));
t=datenum(expdate);


%CREATE DATA FOR BACKGROUND ELECTRIC FIELDS
Exit=zeros(llon,llat,lt);
Eyit=zeros(llon,llat,lt);
for it=1:lt
  Exit(:,:,it)=zeros(llon,llat);   %V/m
  Eyit(:,:,it)=zeros(llon,llat);
end



%CREATE DATA FOR BOUNDARY CONDITIONS FOR POTENTIAL SOLUTION
flagdirich=1;   %if 0 data is interpreted as FAC, else we interpret it as potential
Vminx1it=zeros(llon,llat,lt);
Vmaxx1it=zeros(llon,llat,lt);
Vminx2ist=zeros(llat,lt);
Vmaxx2ist=zeros(llat,lt);
Vminx3ist=zeros(llon,lt);
Vmaxx3ist=zeros(llon,lt);
Etarg=50e-3;            % target E value in V/m
pk=Etarg*sigx2.*xg.h2(lx1,floor(lx2/2),1).*sqrt(pi)./2;
x2ctr=1/2*(xg.x2(lx2)+xg.x2(1));
for it=1:lt
    Vminx1it(:,:,it)=zeros(llon,llat);
    Vmaxx1it(:,:,it)=pk.*erf((MLON-mlonmean)/mlonsig);%.*erf((MLAT-mlatmean)/mlatsig);
     Vminx2ist(:,it)=zeros(llat,1);     %these are just slices
     Vmaxx2ist(:,it)=zeros(llat,1);
     Vminx3ist(:,it)=zeros(llon,1);
     Vmaxx3ist(:,it)=zeros(llon,1);
end



%SAVE THESE DATA TO APPROPRIATE FILES - LEAVE THE SPATIAL AND TEMPORAL INTERPOLATION TO THE
%FORTRAN CODE IN CASE DIFFERENT GRIDS NEED TO BE TRIED.  THE EFIELD DATA DO
%NOT TYPICALLY NEED TO BE SMOOTHED.
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
    filename=[outdir, datelab(ymd,UTsec), '.dat']
    fid=fopen(filename,'w');

    %FOR EACH FRAME WRITE A BC TYPE AND THEN OUTPUT BACKGROUND AND BCs
    fwrite(fid,flagdirich,'real*8');
    fwrite(fid,Exit(:,:,it),'real*8');
    fwrite(fid,Eyit(:,:,it),'real*8');
    fwrite(fid,Vminx1it(:,:,it),'real*8');
    fwrite(fid,Vmaxx1it(:,:,it),'real*8');
    fwrite(fid,Vminx2ist(:,it),'real*8');
    fwrite(fid,Vmaxx2ist(:,it),'real*8');
    fwrite(fid,Vminx3ist(:,it),'real*8');
    fwrite(fid,Vmaxx3ist(:,it),'real*8');

    fclose(fid);
end


%ALSO CREATE A MATLAB OUTPUT FILE FOR GOOD MEASURE
save([outdir,'fields.mat'],'mlon','mlat','MLAT','MLON','Exit','Eyit','Vminx*','Vmax*','expdate');
