%SIMULATIONS LOCAITONS
simname='tohoku20113D_highres_long/';
simname_control='tohoku20113D_highres_long_control/';
basedir='~/zettergmdata/simulations/';
direc=[basedir,simname];
system(['mkdir ',direc,'/TECplots']);    %store output plots with the simulation data


%PATH TO UTILITIES
addpath ../script_utils;


%READ IN THE SIMULATION INFORMATION
[ymd0,UTsec0,tdur,dtout,flagoutput,mloc]=readconfig([direc,'/inputs/config.dat']);


%WE ALSO NEED TO LOAD THE GRID FILE (UNLESS IT ALREADY EXISTS IN THE WORKSPACE)
if (~exist('xg','var'))
  fprintf('Reading grid...\n');
  xg=readgrid([direc,'/inputs/']);
  lx1=xg.lx(1); lx2=xg.lx(2); lx3=xg.lx(3);
  lh=lx1;   %possibly obviated in this version - need to check
  if (lx3==1)
    flag2D=1
  else
    flag2D=0;
  end
end
fprintf('Grid loaded...\n');


%DEFINE A CENTER AND REGION OF INTEREST
if (isempty(mloc))    %in case this run didn't have a disturbance!
  mlatsrc=(pi/2-mean(xg.theta(:)))*180/pi;
  mlonsrc=mean(xg.phi(:))*180/pi;
else
  mlatsrc=mloc(1);
  mlonsrc=mloc(2);
end
thdist=pi/2-mlatsrc*pi/180;    %zenith angle of source location
phidist=mlonsrc*pi/180;


%ANGULAR RANGE TO COVER FOR TEC CALCULATIONS
dang=3.5;
%dang=90;


%TIMES OF INTEREST (MEASURED IN SECONDS FROM BEGINNING SIMJULATION DAY START
times=UTsec0:dtout:UTsec0+tdur;


%NEW (PLOT) GRID SIZE IN R,TH
Re=6370e3;
%lth=250;
lth=1500;
%lr=250;
lr=500;
lphi=250;


%DEFINE A GRID FOR THE INTERPOLATION
rvals=xg.r(1:lh,:,:);
thvals=xg.theta(1:lh,:,:);
phivals=xg.phi(1:lh,:,:);
rmin=min(rvals(:));
rmax=max(rvals(:));
thmin=min(thvals(:));
thmax=max(thvals(:));
phimin=min(phivals(:));
phimax=max(phivals(:));

theta=linspace(thmin,thmax,lth);
r=linspace(rmin,rmax,lr)';
if (~flag2D)       %in case a 2D run has been done, we assume that it is done in x1,x2 (not x1,x3, which would require some changes)
  phi=linspace(phimin,phimax,lphi);
else
  phi=phidist;
end
[THETA,R,PHI]=meshgrid(theta,r,phi);


%THESE ARE DIPOLE COORDINATES OF 
qI=(Re./R).^2.*cos(THETA);
pI=R./Re./sin(THETA).^2;
X3I=PHI;   %phi variable name already used, this is a bit kludgey


ith1=min(find(theta-(thdist-dang*pi/180)>0))
if (isempty(ith1))
   ith1=1; 
end
ith2=min(find(theta-(thdist+dang*pi/180)>0))
if (isempty(ith2))
   ith2=numel(theta); 
end
if (~flag2D)
  iphi1=min(find(phi-(phidist-dang*pi/180)>0))
  if (isempty(iphi1))
    iphi1=1; 
  end
  iphi2=min(find(phi-(phidist+dang*pi/180)>0))
  if (isempty(iphi2))
    iphi2=numel(phi);
  end
else
  iphi1=1;
  iphi2=1;
end
mlat=fliplr(90-theta(ith1:ith2)*180/pi);
mlong=phi(iphi1:iphi2)*180/pi;
itop=lr-1;


%MAIN LOOP FOR TEC CALCULATION
fprintf('Processing %d files...\n',numel(times))
ymd=ymd0;
UTsec=UTsec0;
vTEC=[];
vTEC_control=[];
dvTEC=[];
simdate_series=[];
autoload=0;
for it=1:length(times)
    %LOAD DIST. FILE
    direc=[basedir,simname];
%    filestr=datelab(ymd,UTsec)
%    if (it ~= 1)      %tack on the decimal part
%      filename=[filestr,'.000000.dat']
%    else
%      filename=[filestr,'.000001.dat']
%    end
%    if (flagoutput==1)
%      loadframe3Dcurv;      %note that p from this script is ion composition
%    elseif (flagoutput==2)
%      loadframe3Dcurvavg;
%    else
%      error('Bad output option...')
%    end
%    loadframe_wrapper;
    [ne] = loadframe(direc,UTsec,ymd,UTsec0,ymd0,autoload,flagoutput,mloc,xg);
    simdate=[ymd,UTsec/3600,0,0];    %create a datevec for matlab


    %DEFINE A MESHGRID BASED ON SIMULATION OUTPUT AND DO INTERPOLATION
    if (~flag2D)
      fprintf('3D interpolation...\n')
      x1=xg.x1(3:end-2);
      x2=xg.x2(3:end-2);
      x3=xg.x3(3:end-2);
      [X2,X1,X3]=meshgrid(x2(:),x1(1:lh)',x3(:));   %loadframe overwrites this (sloppy!) so redefine eeach time step
  
      neI=interp3(X2,X1,X3,ne,pI(:),qI(:),X3I(:));
    else
      fprintf('2D interpolation...\n')
      x1=xg.x1(3:end-2);
      x2=xg.x2(3:end-2);
      x3=xg.x3(3:end-2);
      [X2,X1]=meshgrid(x2(:),x1(1:lh)');

      neI=interp2(X2,X1,ne,pI(:),qI(:));
    end


    %RESHAPE AND GET RID OF NANS
    neI=reshape(neI,size(R));
    inds=find(isnan(neI));
    neI(inds)=0;


    %LOAD CONTROL SIMULATION
%    direc=['~/simulations/',simname_control]
    direc=[basedir,simname_control];
%    if (flagoutput==1)
%      loadframe3Dcurv;      %note that p from this script is ion composition
%    elseif (flagoutput==2)
%      loadframe3Dcurvavg;
%    else
%      error('Bad output option...')
%    end
%    loadframe_wrapper;
    [ne] = loadframe(direc,UTsec,ymd,UTsec0,ymd0,autoload,flagoutput,mloc,xg);


    %DEFINE A MESHGRID BASED ON CONTROL SIMULATION OUTPUT AND DO INTERPOLATION
    if (~flag2D)
      fprintf('3D interpolation...\n')
      x1=xg.x1(3:end-2);
      x2=xg.x2(3:end-2);
      x3=xg.x3(3:end-2);
      [X2,X1,X3]=meshgrid(x2(:),x1(1:lh)',x3(:));   %loadframe overwrites this (sloppy!) so redefine eeach time step

      neI_control=interp3(X2,X1,X3,ne,pI(:),qI(:),X3I(:));
    else
      fprintf('2D interpolation...\n')
      x1=xg.x1(3:end-2);
      x2=xg.x2(3:end-2);
      x3=xg.x3(3:end-2);
      [X2,X1]=meshgrid(x2(:),x1(1:lh)');

      neI_control=interp2(X2,X1,ne,pI(:),qI(:));
    end


    %RESHAPE AND GET RID OF NANS IN CONTROL SIMULATION
    neI_control=reshape(neI_control,size(R));
    inds=find(isnan(neI_control));
    neI_control(inds)=0;


    %NOW INTEGRATIOPN TO GET TEC
    if (~flag2D)
      fprintf('Integrating in 3D...\n');
      intne=cumtrapz(r,neI);               %the radial dimension is the first of the neI array
      TECrawnow=intne(itop,ith1:ith2,iphi1:iphi2);
      TECrawnow=squeeze(TECrawnow);        %now the arrays are mlat (1st dim), mlon (2nd dim)
      TECrawnow=flipud(TECrawnow)/1e16;    %mlat runs against x2, scale to TECU
      vTEC=cat(3,vTEC,TECrawnow);          %compile a time series array

      intne=cumtrapz(r,neI_control);               %the radial dimension is the first of the neI array
      TECrawnow=intne(itop,ith1:ith2,iphi1:iphi2);
      TECrawnow=squeeze(TECrawnow);        %now the arrays are mlat (1st dim), mlon (2nd dim)
      TECrawnow=flipud(TECrawnow)/1e16;    %mlat runs against x2, scale to TECU
      vTEC_control=cat(3,vTEC_control,TECrawnow);          %compile a time series array

      dvTEC=cat(3,dvTEC,vTEC(:,:,it)-vTEC_control(:,:,it));
    else
      fprintf('Integrating in 2D...\n');
      intne=cumtrapz(r,neI);               %first dim is radial, as with 3D case
      TECrawnow=intne(itop,ith1:ith2,iphi1:iphi2);
      TECrawnow=squeeze(TECrawnow);        %now the arrays are mlat (1st dim)
      TECrawnow=TECrawnow(:);              %force into a column vector since squeeze will make a row
      TECrawnow=flipud(TECrawnow)/1e16;    %mlat runs against x2, scale to TECU
      vTEC=cat(2,vTEC,TECrawnow);          %compile a time series array

      intne=cumtrapz(r,neI_control);               %first dim is radial, as with 3D case
      TECrawnow=intne(itop,ith1:ith2,iphi1:iphi2);
      TECrawnow=squeeze(TECrawnow);        %now the arrays are mlat (1st dim)
      TECrawnow=TECrawnow(:);              %force into a column vector since squeeze will make a row
      TECrawnow=flipud(TECrawnow)/1e16;    %mlat runs against x2, scale to TECU
      vTEC_control=cat(2,vTEC_control,TECrawnow);          %compile a time series array

      dvTEC=cat(2,dvTEC,vTEC(:,it)-vTEC_control(:,it));    %not vTEC now a 2D array
    end


    %CREATE A DATEVEC FOR THIS SIM TIME SERIES
    simdate_series=[simdate_series;simdate];
    [ymd,UTsec]=dateinc(dtout,ymd,UTsec);

    
    %PLOT THE TOTAL ELECTRON CONTENT EACH TIME FRAME IF WE HAAVE DONE A 3D SIMULATION, OTHERWISE WAIT UNTIL THE END OR A SINGLE PLOT
    if (~flag2D)
      fprintf('Printing TEC plot for current time frame...\n');
      direc=[basedir,simname]
      filename=datelab(ymd,UTsec)
      FS=18;
      imagesc(mlong,mlat,dvTEC(:,:,it));
      colormap(parula(256));
      set(gca,'FontSize',FS);
      axis xy;
      axis tight;
      caxis([-0.25,0.25]);
%      caxis([-4,4]);
      c=colorbar
      set(c,'FontSize',FS)
      xlabel(c,'\Delta vTEC (TECU)')
      xlabel('magnetic long. (deg.)')
      ylabel('magnetic lat. (deg.)')
      hold on;
      ax=axis;
      plot(mlonsrc,mlatsrc,'r^','MarkerSize',10,'LineWidth',2);
      hold off;
      titlestring=datestr(datenum(simdate));
      title(titlestring);
      print('-dpng',[direc,'/TECplots/',filename,'.png'],'-r300');
    end
end


%CREATE A FULL TIME SERIES PLOT IF IN 2D
if (flag2D)
  fprintf('Printing TEC plot for entire time series...\n');
  direc=[basedir,simname]
  FS=18;
  t=datenum(simdate_series);
  imagesc(t,mlat,dvTEC(:,:));
  colormap(parula(256));
  set(gca,'FontSize',FS);
  axis xy;
  datetick;
  axis tight;
  caxis([-0.25,0.25]);
  c=colorbar
  set(c,'FontSize',FS)
  xlabel(c,'\Delta vTEC (TECU)')
  xlabel('UT')
  ylabel('magnetic lat. (deg.)')
  hold on;
  ax=axis;
  plot([t(1), t(end)],[mlatsrc, mlatsrc],'r--','MarkerSize',10,'LineWidth',2);
  hold off;
  print('-dpng',[direc,'/TECplots/TEC_timeseries.png'],'-r300');
end


%SAVE THE DATA TO A .MAT FILE IN CASE WE3 NEED IT LATER
t=datenum(simdate_series);
save([direc,'/vTEC.mat'],'mlat','mlong','t','simdate_series','*vTEC*','-v7');


%RESET Path
rmpath ../script_utils;
