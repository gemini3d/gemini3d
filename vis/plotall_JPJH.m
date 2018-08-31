%SIMULATIONS LOCAITONS
%simname='curvtest_tohoku_medres_2D_electro/';
%simname='curvtest_tohoku_highres_weak/';
%simname='curvtest_tohoku_highres_2D/'
%simname='curvtest_tohoku_medres_2D_geomag/'
%simname='curvtest_tohoku_medres_2D_neuinfo/'
%simname='chile2015/'
%simname='mooreOK/'
%simname='3DPCarc/'
%simname='GDI_periodic_highres/'
%simname='KHI_periodic_highres_fileinput/'
%simname='chile20153D/'
%simname='GDI_nonperiodic/'
%simname='isinglass/'
%simname='chile2015_0.5_strong/'
%simname='nepal2015/'
%simname='chile20153D_0.5_medhighres/'
%simname='2Dtest_fields/'
simname='ARCS_current2/'


%MAKE DIRECTORIES FOR OUTPUT FILES
%basedir='/Volumes/TravelDisk2/simulations/'
basedir='~/zettergmdata/simulations/'
direc=[basedir,simname];
system(['mkdir ',direc,'/nplots']);    %store output plots with the simulation data
system(['mkdir ',direc,'/v1plots']);    %store output plots with the simulation data
system(['mkdir ',direc,'/v2plots']);    %store output plots with the simulation data
system(['mkdir ',direc,'/v3plots']);    %store output plots with the simulation data
system(['mkdir ',direc,'/J1plots']);    %store output plots with the simulation data
system(['mkdir ',direc,'/Tiplots']);    %store output plots with the simulation data
system(['mkdir ',direc,'/Teplots']);    %store output plots with the simulation data
system(['mkdir ',direc,'/J2plots']);    %store output plots with the simulation data
system(['mkdir ',direc,'/J3plots']);    %store output plots with the simulation data
system(['mkdir ',direc,'/JPplots']);
system(['mkdir ',direc,'/JHplots']);


%PATH TO PLOTTIN FUNCTIONS AND SHARED SCRIPT UTILITIES
addpath ./plotfunctions;
addpath ../script_utils;


%%READ IN THE SIMULATION INFORMATION
[ymd0,UTsec0,tdur,dtout,flagoutput,mloc]=readconfig([direc,'/inputs/config.dat']);



%%CHECK WHETHER WE NEED TO RELOAD THE GRID (WHICH CAN BE TIME CONSUMING)
%if (~exist('xg','var'))
%  %WE ALSO NEED TO LOAD THE GRID FILE
%  xg=readgrid([direc,'/']);
%  fprintf('Grid loaded...\n');
%end
%

%FUNCTINO TO PLOT WITH
%plotfun=@plot2D_curv;
%plotfun=@plot2D_cart;
%plotfun=@plot2D_curv_north;
%plotfun=@plot2D_curv_south;
%plotfun=@plot3D_curv_frames;
%plotfun=@plot3D_cart_frames;
%plotfun=@plot3D_curv_frames_long;
plotfun=@plot3D_cart_frames_long_ENU;


%%COMPUTE SOURUCE LOCATION IN MCOORDS
%if (~isempty(mloc))
%  mlat=mloc(1);
%  mlon=mloc(2);
%else
%  mlat=[];
%  mlon=[];
%end



%TIMES OF INTEREST
times=UTsec0:dtout:UTsec0+tdur;


ymd=ymd0;
UTsec=UTsec0;
for it=1:length(times)
%    %LOAD DIST. FILE
%    filestr=datelab(ymd,UTsec)
%    if (it ~= 1)      %tack on the decimal part
%      filename=[filestr,'.000000.dat']
%    else
%      filename=[filestr,'.000001.dat']
%    end
%%    if(xg.lx(2)==1 | xg.lx(3)==1)
%%      loadframe3Dcurv;
%%    else
%%      loadframe3Dcurvavg;
%%    end
%    if (flagoutput==1)
%      loadframe3Dcurv;
%    elseif (flagoutput==2)
%      loadframe3Dcurvavg;
%    else
%      loadframe3Dcurvne;
%    end
  loadframe_wrapper;

%
%  %RESOLVE CURRENTS INTO PEDERSEN AND HALL COMPONENTS
%  lx1=xg.lx(1); lx2=xg.lx(2); lx3=xg.lx(3);
%  vperp2=permute(squeeze(v2(floor(lx1/2),:,:)),[2,1]);
%%  vperp2=reshape(vperp2,[lx3,lx2,1]);  
%  vperp2=repmat(vperp2,[1,1,lx1]);
%  vperp3=permute(squeeze(v3(floor(lx1/2),:,:)),[2,1]);
%%  vperp3=reshape(vperp3,[lx3,lx2,1]);
%  vperp3=repmat(vperp3,[1,1,lx1]);
%  v=cat(4,vperp3,vperp2,zeros(lx3,lx2,lx1));    %ExB drift
%  magvperp=sqrt(sum(v.^2,4));
%  evperp=v./repmat(magvperp,[1,1,1,3]);         %unit vector for ExB drift
%  e1curv=cat(4,zeros(lx3,lx2,lx1),zeros(lx3,lx2,lx1),ones(lx3,lx2,lx1));   %unit vector in curvilinear basis
%  B=cat(4,zeros(lx3,lx2,lx1),zeros(lx3,lx2,lx1),ones(lx3,lx2,lx1).*permute(xg.Bmag,[3,2,1]));
%  E=cross(-1*v,B,4);        %electric field
%  magE=sqrt(sum(E.^2,4));
%  eE=E./repmat(magE,[1,1,1,3]);
%
%  J=cat(4,permute(J3,[3,2,1]),permute(J2,[3,2,1]),permute(J1,[3,2,1]));    %current density vector
%  JH=dot(J,-1*evperp,4);
%  JP=dot(J,eE,4);
%  Jfac=cat(4,zeros(lx3,lx2,lx1),zeros(lx3,lx2,lx1),permute(J1,[3,2,1]));
%

  %RESOLVE CURRENTS INTO PEDERSEN AND HALL COMPONENTS
  lx1=xg.lx(1); lx2=xg.lx(2); lx3=xg.lx(3);
  vperp2=squeeze(v2(floor(lx1/2),:,:));
  vperp2=repmat(vperp2,[1,1,lx1]);
  vperp3=squeeze(v3(floor(lx1/2),:,:));
  vperp3=repmat(vperp3,[1,1,lx1]);
  v=cat(4,vperp2,vperp3,zeros(lx2,lx3,lx1));    %ExB drift
  magvperp=sqrt(sum(v.^2,4));
  evperp=v./repmat(magvperp,[1,1,1,3]);         %unit vector for ExB drift
  e1curv=cat(4,zeros(lx2,lx3,lx1),zeros(lx2,lx3,lx1),ones(lx2,lx3,lx1));   %unit vector in curvilinear basis
  B=cat(4,zeros(lx2,lx3,lx1),zeros(lx2,lx3,lx1),ones(lx2,lx3,lx1).*permute(xg.Bmag,[2,3,1]));
  E=cross(-1*v,B,4);        %electric field
  magE=sqrt(sum(E.^2,4));
  eE=E./repmat(magE,[1,1,1,3]);

  J=cat(4,permute(J2,[2,3,1]),permute(J3,[2,3,1]),permute(J1,[2,3,1]));    %current density vector
  JH=dot(J,-1*evperp,4);     %projection of current in -evperp direction
  JHvec=-1*evperp.*repmat(JH,[1,1,1,3]);
  JP=dot(J,eE,4);    %project of current in eE unit vector direction
  JPvec=eE.*repmat(JP,[1,1,1,3]);
  Jfac=cat(4,zeros(lx2,lx3,lx1),zeros(lx2,lx3,lx1),permute(J1,[2,3,1]));    %field aligned current vector

  %MAKE PLOTS OF OUR FAVORITE PARAMETERS
  figure(10);
  clf;
%  h2=plotfun(ymd,UTsec,xg,ne(:,:,:),'n_e (m^{-3})',[0 8e11],[mlatsrc,mlonsrc]);
  h2=plotfun(ymd,UTsec,xg,ne(:,:,:),'n_e (m^{-3})',[0 2e11],[mlatsrc,mlonsrc]);
  print('-dpng',[direc,'/nplots/',filename,'.png'],'-r300');

  if (flagoutput~=3)
    figure(1);
    clf;
    h2=plotfun(ymd,UTsec,xg,v1(:,:,:),'v_1 (m/s)',[-75 75],[mlatsrc,mlonsrc]);
    print('-dpng',[direc,'/v1plots/',filename,'.png'],'-r300');
   
    figure(2);
    clf;
    h2=plotfun(ymd,UTsec,xg,Ti(:,:,:),'T_i (K)',[100 2500],[mlatsrc,mlonsrc]);
    print('-dpng',[direc,'/Tiplots/',filename,'.png'],'-r300');
  
    figure(3);
    clf;
    h2=plotfun(ymd,UTsec,xg,Te(:,:,:),'T_e (K)',[100 6000],[mlatsrc,mlonsrc]);
    print('-dpng',[direc,'/Teplots/',filename,'.png'],'-r300');
   
    figure(4);
    clf;
%    h2=plotfun(ymd,UTsec,xg,J1(:,:,:)*1e6,'J_1 (uA/m^2)',[-0.025 0.025],[mlatsrc,mlonsrc]);
    h2=plotfun(ymd,UTsec,xg,J1(:,:,:)*1e6,'J_1 (uA/m^2)',[-10 10],[mlatsrc,mlonsrc]);
    print('-dpng',[direc,'/J1plots/',filename,'.png'],'-r300'); 
   
    figure(5);
    clf;
    h2=plotfun(ymd,UTsec,xg,v2(:,:,:),'v_2 (m/s)',[-500 500],[mlatsrc,mlonsrc]);
    print('-dpng',[direc,'/v2plots/',filename,'.png'],'-r300'); 
   
    figure(6);
    clf;
    h2=plotfun(ymd,UTsec,xg,v3(:,:,:),'v_3 (m/s)',[-500 500],[mlatsrc,mlonsrc]);
    print('-dpng',[direc,'/v3plots/',filename,'.png'],'-r300');
    
    figure(7);
    clf;
    h2=plotfun(ymd,UTsec,xg,J2(:,:,:)*1e6,'J_2 (uA/m^2)',[-10 10],[mlatsrc,mlonsrc]);
    print('-dpng',[direc,'/J2plots/',filename,'.png'],'-r300'); 
   
    figure(8);
    clf;
    h2=plotfun(ymd,UTsec,xg,J3(:,:,:)*1e6,'J_3 (uA/m^2)',[-10 10],[mlatsrc,mlonsrc]);
    print('-dpng',[direc,'/J3plots/',filename,'.png'],'-r300');    

    figure(11);
    clf;
    JPplot=sqrt(sum(JPvec.^2,4));
    JPplot=ipermute(JPplot,[2,3,1]);
    h2=plotfun(ymd,UTsec,xg,log10(JPplot(:,:,:)),'log_{10} |J_P| (A/m^2)',[-7 -5],[mlatsrc,mlonsrc]);
    print('-dpng',[direc,'/JPplots/',filename,'.png'],'-r300');

    figure(12);
    clf;
    JHplot=sqrt(sum(JHvec.^2,4));
    JHplot=ipermute(JHplot,[2,3,1]);
    h2=plotfun(ymd,UTsec,xg,log10(JHplot(:,:,:)),'|J_H| (log_{10} A/m^2)',[-7 -5],[mlatsrc,mlonsrc]);
    print('-dpng',[direc,'/JHplots/',filename,'.png'],'-r300');
  end

  %close all;
  [ymd,UTsec]=dateinc(dtout,ymd,UTsec);
end


%MAKE A BINARY FILE WITH THE FINAL TIME STEP OUTPUT OF HALL AND PEDERSEN CURRENTS
%JP=permute(JP,[2,1,3,4]);
%JH=permute(JH,[2,1,3,4]);
%Jfac=permute(Jfac,[2,1,3,4]);
fid=fopen('JPJH.dat','w');
datsize=[lx2,lx3,lx1];
fwrite(fid,datsize,'integer*4');
fwrite(fid,xg.x2(3:end-2),'real*8');
fwrite(fid,xg.x3(3:end-2),'real*8');
fwrite(fid,xg.x1(3:end-2),'real*8');
fwrite(fid,JPvec,'real*8');
fwrite(fid,JHvec,'real*8');
fwrite(fid,Jfac,'real*8');
fclose(fid);


%CLEAR THE PATHS ADDED
rmpath ./plotfunctions;
rmpath ../script_utils;

