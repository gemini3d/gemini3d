%SIMULATIONS LOCAITONS
simname='curvtest_tohoku_medres_2D_electro/';


%MAKE DIRECTORIES FOR OUTPUUT FILES
direc=['~/simulations/',simname];
system(['mkdir ',direc,'/nplots']);    %store output plots with the simulation data
system(['mkdir ',direc,'/v1plots']);    %store output plots with the simulation data
system(['mkdir ',direc,'/v2plots']);    %store output plots with the simulation data
system(['mkdir ',direc,'/v3plots']);    %store output plots with the simulation data
system(['mkdir ',direc,'/J1plots']);    %store output plots with the simulation data
system(['mkdir ',direc,'/Tiplots']);    %store output plots with the simulation data
system(['mkdir ',direc,'/Teplots']);    %store output plots with the simulation data


%DATE AND FILE STRINGS
%ymd0=[2011,03,10];
%UTsec0=20783;
ymd0=[2011,03,11];
UTsec0=20783;


%CHECK WHETHER WE NEED TO RELOAD THE GRID
if (~exist('xg','var'))
  %WE ALSO NEED TO LOAD THE GRID FILE
  xg=readgrid([direc,'/']);
  fprintf('Grid loaded...\n');
end


%FUNCTINO TO PLOT WITH
plotfun=@plot2D_curv;


%COMPUTE SOURUCE LOCATIOKN IN MCOORDS
glat=38.429575d0
glon=142.734757d0
addpath ../setup;
[theta,phi]=geog2geomag(glat,glon);
mlat=90-theta*180/pi
mlon=phi*180/pi
rmpath ../setup;


%TIMES OF INTEREST
%dtout=1800;
%times=UTsec0:dtout:UTsec0+86400;
dtout=30;
times=UTsec0:dtout:UTsec0+3600;


ymd=ymd0;
UTsec=UTsec0;
for it=1:length(times)
    %LOAD DIST. FILE
    filestr=datelab(ymd,UTsec)
    if (it ~= 1)      %tack on the decimal part
      filename=[filestr,'.000000.dat']
    else
      filename=[filestr,'.000001.dat']
    end
    loadframe3Dcurv;


  %MAKE PLOTS OF OUR FAVORITE PARAMETERS
  figure;
  h2=plotfun(ymd,UTsec,xg,log10(ns(:,:,:,7)),'log_{10} n_e (m^{-3})',[8 12.5],[mlat,mlon]);
  print('-dpng',[direc,'/nplots/',filename,'.png'],'-r300');

  figure;
  h2=plotfun(ymd,UTsec,xg,v1(:,:,:),'v_1 (m/s)',[-150 150],[mlat,mlon]);
  print('-dpng',[direc,'/v1plots/',filename,'.png'],'-r300');
 
  figure;
  h2=plotfun(ymd,UTsec,xg,Ti(:,:,:),'T_i (K)',[100 2500],[mlat,mlon]);
  print('-dpng',[direc,'/Tiplots/',filename,'.png'],'-r300');

  figure;
  h2=plotfun(ymd,UTsec,xg,Ts(:,:,:,7),'T_e (K)',[100 6000],[mlat,mlon]);
  print('-dpng',[direc,'/Teplots/',filename,'.png'],'-r300');
 
  figure;
  h2=plotfun(ymd,UTsec,xg,J1(:,:,:)*1e6,'J_1 (uA/m^2)',[-0.025 0.025],[mlat,mlon]);
  print('-dpng',[direc,'/J1plots/',filename,'.png'],'-r300'); 
 
  figure;
  h2=plotfun(ymd,UTsec,xg,v2(:,:,:),'v_2 (m/s)',[-5 5],[mlat,mlon]);
  print('-dpng',[direc,'/v2plots/',filename,'.png'],'-r300'); 
 
  figure;
  h2=plotfun(ymd,UTsec,xg,v3(:,:,:),'v_3 (m/s)',[-5 5],[mlat,mlon]);
  print('-dpng',[direc,'/v3plots/',filename,'.png'],'-r300');

  close all;
  [ymd,UTsec]=dateinc(dtout,ymd,UTsec);
end
