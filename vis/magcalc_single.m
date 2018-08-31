%SIMULATIONS LOCAITONS
simname='chile20153D/';
%direc=['~/simulations/',simname];
direc=['/Volumes/TravelDisk2/simulations/',simname];
system(['mkdir ',direc,'/magplots']);    %store output plots with the simulation data


%UTseconds of the frame of interest
ymd_TOI=[2015,09,16];
UTsec_TOI=82923;


%ADD PATHS
addpath ../script_utils;


%SIMULATION META-DATA
[ymd0,UTsec0,tdur,dtout,flagoutput,mloc]=readconfig([direc,'/inputs/config.dat']);


%TABULATE THE SOURCE LOCATION
mlatsrc=mloc(1);
mlonsrc=mloc(2);
thdist=pi/2-mlatsrc*pi/180;    %zenith angle of source location
phidist=mlonsrc*pi/180;


%ANGULAR RANGE TO COVER FOR THE CALCLUATIONS (THIS IS FOR THE FIELD POINTS - SOURCE POINTS COVER ENTIRE GRID)
dang=5;


%WE ALSO NEED TO LOAD THE GRID FILE
if (~exist('xg','var'))
  fprintf('Reading grid...\n');
  xg=readgrid([direc,'/']);
  lx1=xg.lx(1); lx2=xg.lx(2); lx3=xg.lx(3);
  lh=lx1;   %possibly obviated in this version - need to check
  if (lx3==1)
    flag2D=1;
  else
    flag2D=0;
  end
end
fprintf('Grid loaded...\n');


%TIMES OF INTEREST (MEASURED IN SECONDS FROM BEGINNING SIMJULATION DAY START
times=UTsec_TOI;
lt=numel(times);



%NEW SOURCE GRID SIZE
Re=6370e3;
lthetap=500;
lrp=250;
%lthetap=200;
%lrp=200;
if (~flag2D)
%  lphip=250;
  lphip=150;
else
  lphip=1;
end

%DEFINE A GRID FOR THE INTERPOLATION ONTO SOURCE COORDINATES (PRIMED COORDS.)
rvals=xg.r(1:lh,:,:);
thvals=xg.theta(1:lh,:,:);
phivals=xg.phi(1:lh,:,:);
rmin=min(rvals(:));
rmax=max(rvals(:));
thmin=min(thvals(:));
thmax=max(thvals(:));
phimin=min(phivals(:));
phimax=max(phivals(:));

thetap=linspace(thmin,thmax,lthetap);
rp=linspace(rmin,rmax,lrp)';
if (~flag2D)       %in case a 2D run has been done, we assume that it is done in x1,x2 (not x1,x3, which would require some changes)
  phip=linspace(phimin,phimax,lphip);
else
  phip=phidist;     %in 2D just assume that we are in the meridian of disturbance
end
[THETAP,RP,PHIP]=meshgrid(thetap,rp,phip);


%THESE ARE DIPOLE COORDINATES OF THE SOURCE POINTS (NEEDED TO PERFORM INTERPOLATION)
QP=(Re./RP).^2.*cos(THETAP);
PP=RP./Re./sin(THETAP).^2;
%phi already handled above - i.e. same as in spherical coordinates
%PHIP;   %phi variable name already used


%FIELD POINTS OF INTEREST (CAN/SHOULD BE DEFINED INDEPENDENT OF SIMULATION GRID)
ltheta=150;
if (~flag2D)
  lphi=150;
else
  lphi=1;
end
%ltheta=100;
%lphi=100;
lr=1;

thmin=thdist-dang*pi/180;
thmax=thdist+dang*pi/180;
phimin=phidist-dang*pi/180;
phimax=phidist+dang*pi/180;

theta=linspace(thmin,thmax,ltheta);
if (~flag2D)
  phi=linspace(phimin,phimax,lphi);
else
  phi=phidist;
end
r=6370e3;                          %use ground level for altitude


%CONVEY THE SOURCE AND FIELD POINTS INTO ACTUAL DISTANCES SO WE CAN ASSUME A CARTESIAN INTEGRATION (THIS IS APPROXIMATE)
rref=mean(r(:));
thetaref=mean(theta(:));
x=r;
y=rref.*theta;
z=rref.*sin(thetaref).*phi;
lx=numel(x);
ly=numel(y);
lz=numel(z);


%MAGNETIC LATITUDE AND LONGITUDE OF FIELD POINTS FOR PLOTTING
mlat=90-theta*180/pi;
mlon=phi*180/pi;
[mlat,imlat]=sort(mlat);
[mlon,imlon]=sort(mlon);   %note that magnetic field data will need to be resorted following calcualtions


%"CARTESIAN" COORDINATES OF SOURCE POINTS (CALCULATED BASED ON FIELD POINT SPECIFICATIONS)
rpref=rref;    %use same reference value as the field coordinates
thetapref=thetaref;
xp=rp;
yp=(rpref).*thetap;
zp=(rpref).*sin(thetapref).*phip;
lxp=numel(xp);
lyp=numel(yp);
lzp=numel(zp);

[YP,XP,ZP]=meshgrid(yp,xp,zp);   %to be consistent with RP,THETAP,PHIP DEFINED ABOVE

dxp=xp(2)-xp(1);     %assume uniform spacing
dyp=yp(2)-yp(1);
if (~flag2D)         %in 2D we don't have this third dimension
  dzp=zp(2)-zp(1);
end


%MAIN TIME LOOP
fprintf('Started processing %d files...\n',numel(times))
ymd=ymd_TOI;
UTsec=UTsec_TOI;
simdate_series=[];
Brt=zeros(lx,ly,lz,lt);
Bthetat=zeros(lx,ly,lz,lt);
Bphit=zeros(lx,ly,lz,lt);
for it=1:lt
    %LOAD DIST. FILE
    filestr=datelab(ymd,UTsec)
%    if (it ~= 1)      %tack on the decimal part
      filename=[filestr,'.000000.dat']
%    else
%      filename=[filestr,'.000001.dat']
%    end
    if (flagoutput==1)
      loadframe3Dcurv
    elseif (flagoutput==2)
      loadframe3Dcurvavg;
    else
      error('Wrong output option...');
    end

    %DEFINE A MESHGRID BASED ON SIMULATION OUTPUT AND DO INTERPOLATION
    if (~flag2D)
      fprintf('3D meshgrid...\n')
      x1=xg.x1(3:end-2);
      x2=xg.x2(3:end-2);
      x3=xg.x3(3:end-2);
      [X2,X1,X3]=meshgrid(x2(:),x1(1:lh)',x3(:));   %loadframe overwrites this (sloppy!) so redefine eeach time step
    else
      fprintf('2D meshgrid...\n')
      x1=xg.x1(3:end-2);
      x2=xg.x2(3:end-2);
      x3=xg.x3(3:end-2);
      [X2,X1]=meshgrid(x2(:),x1(1:lh)');
    end


    %MAKE SURE THERE ARE NO CURRENTS BELOW THE IONOSPHERE
    inds=find(xg.alt<80e3);
    J1(inds)=0;  J2(inds)=0; J3(inds)=0;


    %ROTATE MODEL DATA INTO SPHERICAL COORDINATES
    fprintf('2D/3D rotations...\n')
    Jr=J1.*dot(xg.er,xg.e1,4)+J2.*dot(xg.er,xg.e2,4)+J3.*dot(xg.er,xg.e3,4);
    Jtheta=J1.*dot(xg.etheta,xg.e1,4)+J2.*dot(xg.etheta,xg.e2,4)+J3.*dot(xg.etheta,xg.e3,4);
    Jphi=J1.*dot(xg.ephi,xg.e1,4)+J2.*dot(xg.ephi,xg.e2,4)+J3.*dot(xg.ephi,xg.e3,4);


    %INTEROPLATE ONTO SOURCE COORDINATES AS PREP FOR INTEGRATIONS
    if (~flag2D)
      fprintf('3D interpolation onto source coords....\n')
      JrI=interp3(X2,X1,X3,Jr,PP(:),QP(:),PHIP(:));
      JthetaI=interp3(X2,X1,X3,Jtheta,PP(:),QP(:),PHIP(:));
      JphiI=interp3(X2,X1,X3,Jphi,PP(:),QP(:),PHIP(:));
    else
      fprintf('2D interpolation onto source coords....\n')
      JrI=interp2(X2,X1,Jr,PP(:),QP(:));
      JthetaI=interp2(X2,X1,Jtheta,PP(:),QP(:));
      JphiI=interp2(X2,X1,Jphi,PP(:),QP(:));
    end
    JrI=reshape(JrI,size(PP));
    inds=find(isnan(JrI));
    JrI(inds)=0;

    JthetaI=reshape(JthetaI,size(PP));
    inds=find(isnan(JthetaI));
    JthetaI(inds)=0;

    JphiI=reshape(JphiI,size(PP));
    inds=find(isnan(JphiI));
    JphiI(inds)=0;


    %RE-ENVISION THE SOURCES AND SOURCE COORDINATES AS LOCAL CARTESIAN xprime,yprime,zprime COORDIANTES
    JxI=JrI;
    JyI=JthetaI;
    JzI=JphiI;

%{
    JxI=JrI;
    JyI=JthetaI;
    JzI=JphiI;
    XP=RP;
    YP=RP.*THETAP;
    ZP=RP.*sin(THETAP).*PHIP;
    %the above coordiantes are not plaid - do we just need to do a proper integration in spherical or convert to Cartesian???  This would require an extra rotations back into geomagnetic once the calculations is done, but it would be "proper" (except for account for the Earth's boundary surface)...
%}


    %PRINT SOME SANITY CHECK DIAGNOSTICS
    min(JxI(:)),max(JxI(:))
    min(JyI(:)),max(JyI(:))
    min(JzI(:)),max(JzI(:))

    min(J1(:)),max(J1(:))
    min(J2(:)),max(J2(:))
    min(J3(:)),max(J3(:))


    %LOOP OVER ***FIELD POINTS*** AND COMPUTE THE MAGNETIC FIELD FOR EACH ONE
    fprintf('Beginning magnetic field integrations for time step:  %d...\n',it);
    mu0=4*pi*1e-7;    
    Bx=zeros(lx,ly,lz);
    By=zeros(lx,ly,lz);
    Bz=zeros(lx,ly,lz);
    for ix=1:lx
      for iy=1:ly
        for iz=1:lz
          fprintf('Integrating field point (%d,%d,%d) of (%d,%d,%d) for time step %d...\n',ix,iy,iz,lx,ly,lz,it);


          %COMPUTE DISTANCES NEEDED FOR INTEGRALS
          Rx=x(ix)-XP;
          Ry=y(iy)-YP;
          Rz=z(iz)-ZP;
          if (~flag2D)
            Rcubed=(Rx.^2+Ry.^2+Rz.^2).^(3/2);   %may need a floor here if the field points include source locations, otherwise will get a singularity
          else
            denom=Rx.^2+Ry.^2;
          end

          %BX COMPONENT
          if (~flag2D)
            intgrnd=mu0/4/pi*(JyI.*Rz-JzI.*Ry)./Rcubed;
            intgrndavg=1/8*( intgrnd(1:lrp-1,1:lthetap-1,1:lphip-1) + intgrnd(2:lrp,1:lthetap-1,1:lphip-1) + ...
                           intgrnd(1:lrp-1,2:lthetap,1:lphip-1) + intgrnd(2:lrp,2:lthetap,1:lphip-1) + ...
                           intgrnd(1:lrp-1,1:lthetap-1,2:lphip) + intgrnd(2:lrp,1:lthetap-1,2:lphip) + ...
                           intgrnd(1:lrp-1,2:lthetap,2:lphip) + intgrnd(2:lrp,2:lthetap,2:lphip) );     %the integrandd averaged to the source coordiante cell center
            Bx(ix,iy,iz)=sum(intgrndavg(:)*dxp*dyp*dzp);
          else
            intgrnd=mu0/4/pi*(-2*JzI.*Ry)./denom;
            intgrndavg=1/4*( intgrnd(1:lrp-1,1:lthetap-1) + intgrnd(2:lrp,1:lthetap-1) + ...
                           intgrnd(1:lrp-1,2:lthetap) + intgrnd(2:lrp,2:lthetap) );
            Bx(ix,iy,iz)=sum(intgrndavg(:)*dxp*dyp);
          end


          %BY
          if (~flag2D)
            intgrnd=-1*mu0/4/pi*(JxI.*Rz-JzI.*Rx)./Rcubed;
            intgrndavg=1/8*( intgrnd(1:lrp-1,1:lthetap-1,1:lphip-1) + intgrnd(2:lrp,1:lthetap-1,1:lphip-1) + ...
                           intgrnd(1:lrp-1,2:lthetap,1:lphip-1) + intgrnd(2:lrp,2:lthetap,1:lphip-1) + ...
                           intgrnd(1:lrp-1,1:lthetap-1,2:lphip) + intgrnd(2:lrp,1:lthetap-1,2:lphip) + ...
                           intgrnd(1:lrp-1,2:lthetap,2:lphip) + intgrnd(2:lrp,2:lthetap,2:lphip) );
            By(ix,iy,iz)=sum(intgrndavg(:)*dxp*dyp*dzp);
          else
            intgrnd=mu0/4/pi*(2*JzI.*Rx)./denom;
            intgrndavg=1/4*( intgrnd(1:lrp-1,1:lthetap-1) + intgrnd(2:lrp,1:lthetap-1) + ...
                           intgrnd(1:lrp-1,2:lthetap) + intgrnd(2:lrp,2:lthetap) );
            By(ix,iy,iz)=sum(intgrndavg(:)*dxp*dyp);
          end


          %BZ
          if (~flag2D)
            intgrnd=mu0/4/pi*(JxI.*Ry-JyI.*Rx)./Rcubed;
            intgrndavg=1/8*( intgrnd(1:lrp-1,1:lthetap-1,1:lphip-1) + intgrnd(2:lrp,1:lthetap-1,1:lphip-1) + ...
                           intgrnd(1:lrp-1,2:lthetap,1:lphip-1) + intgrnd(2:lrp,2:lthetap,1:lphip-1) + ...
                           intgrnd(1:lrp-1,1:lthetap-1,2:lphip) + intgrnd(2:lrp,1:lthetap-1,2:lphip) + ...
                           intgrnd(1:lrp-1,2:lthetap,2:lphip) + intgrnd(2:lrp,2:lthetap,2:lphip) );
            Bz(ix,iy,iz)=sum(intgrndavg(:)*dxp*dyp*dzp);
          else
            intgrnd=mu0/4/pi*2*(JxI.*Ry-JyI.*Rx)./denom;
            intgrndavg=1/4*( intgrnd(1:lrp-1,1:lthetap-1) + intgrnd(2:lrp,1:lthetap-1) + ...
                           intgrnd(1:lrp-1,2:lthetap) + intgrnd(2:lrp,2:lthetap) );
            Bz(ix,iy,iz)=sum(intgrndavg(:)*dxp*dyp);
          end
        end
      end
    end


    %SORT THE MAGNETIC FIELD CALCULATIONS ACCORDING OT MLAT AND MLON
    Bx=Bx(:,imlat,:);    %not sure these commute so I'll do things in steps...
    Bx=Bx(:,:,imlon);
    By=By(:,imlat,:);
    By=By(:,:,imlon);
    Bz=Bz(:,imlat,:);
    Bz=Bz(:,:,imlon);


    %RENAME THE MAGNETIC FIELD COMPONENTS ACCORDING TO ACTUALY GEOMETRY USED (SHERICAL) and convert to nanotesla
    Br=Bx;
    Btheta=By;
    Bphi=Bz;
    Brt(:,:,:,it)=Br;
    Bthetat(:,:,:,it)=Btheta;
    Bphit(:,:,:,it)=Bphi;


    %GET RID OF ANY SINGLETON DIMENSIONS (in 3d this reduces the vertical dim, in 2D vertical and longtude)
    Br=squeeze(Br);
    Btheta=squeeze(Btheta);
    Bphi=squeeze(Bphi);


    %CONSTRUCT A RUNNING DATE VECTOR, IN CASE IT IS USEFUL IN LATER CALCULATIONS
    simdate_series=[simdate_series;simdate];
    [ymd,UTsec]=dateinc(dtout,ymd,UTsec);

    
    %PLOT THE TOTAL ELECTRON CONTENT OF THIS FRAEM IF THIS WAS A 3D RUN, OTHERWISE WAIT UNTIL THE END TO DO ONE TIME SERIES PLOT
    if (~flag2D)
      fprintf('Printing magnetic field plots...\n');
      figure;
      set(gcf,'PaperPosition',[0 0 11 4]);
      FS=10;

      titlestring=datestr(datenum(simdate));
  
      subplot(131);
      param=Br*1e9;
      imagesc(mlon,mlat,param);
      colormap(parula(256));
      set(gca,'FontSize',FS);
      axis xy;
      axis tight;
      caxlim=max(abs(param(:)))
%      caxlim=max(caxlim,0.05)
      caxlim=max(caxlim,0.001);
      caxis([-caxlim,caxlim]);
      c=colorbar
      set(c,'FontSize',FS)
  %    xlabel(c,'B_r (nT)')
      title(['B_r (nT)  ',titlestring])
      xlabel('magnetic long. (deg.)')
      ylabel('magnetic lat. (deg.)')
      hold on;
      ax=axis;
      plot(mlonsrc,mlatsrc,'r^','MarkerSize',10,'LineWidth',2);
      hold off;
  
      subplot(132);
      param=Btheta*1e9;
      imagesc(mlon,mlat,param);
      colormap(parula(256));
      set(gca,'FontSize',FS);
      axis xy;
      axis tight;
      caxlim=max(abs(param(:)))
%      caxlim=max(caxlim,0.05)
      caxlim=max(caxlim,0.001);
      caxis([-caxlim,caxlim]);
      c=colorbar
      set(c,'FontSize',FS)
  %    xlabel(c,'B_\theta (nT)')
      title('B_\theta (nT)')
      xlabel('magnetic long. (deg.)')
      ylabel('magnetic lat. (deg.)')
      hold on;
      ax=axis;
      plot(mlonsrc,mlatsrc,'r^','MarkerSize',10,'LineWidth',2);
      hold off;
  
      subplot(133);
      param=Bphi*1e9;
      imagesc(mlon,mlat,param);
      colormap(parula(256));
      set(gca,'FontSize',FS);
      axis xy;
      axis tight;
      caxlim=max(abs(param(:)))
%      caxlim=max(caxlim,0.05)
      caxlim=max(caxlim,0.001);
      caxis([-caxlim,caxlim]);
      c=colorbar
      set(c,'FontSize',FS)
  %    xlabel(c,'B_\phi (nT)')
      title('B_\phi (nT)')
      xlabel('magnetic long. (deg.)')
      ylabel('magnetic lat. (deg.)')
      hold on;
      ax=axis;
      plot(mlonsrc,mlatsrc,'r^','MarkerSize',10,'LineWidth',2);
      hold off;
  
  
      print('-dpng',[direc,'/magplots/',filename,'.png'],'-r300');
      close all;
    end
end


%IF IT WAS A 2D RUN PRODUCE A FULL TIME SERIES PLOT AFTERWARD
if (flag2D)
  Br=squeeze(Brt);    %dimensions of mlat,time after squeeze
  Btheta=squeeze(Bthetat);
  Bphi=squeeze(Bphit);
  t=datenum(simdate_series);

  fprintf('Printing magnetic field plots...\n');
  figure;
  set(gcf,'PaperPosition',[0 0 11 4]);
  FS=10;
  
  subplot(131);
  param=Br*1e9;
  imagesc(t,mlat,param);
  colormap(parula(256));
  set(gca,'FontSize',FS);
  axis xy;
  datetick;
  axis tight;
  caxlim=max(abs(param(:)))
%  caxlim=max(caxlim,0.05)
  caxis([-caxlim,caxlim]);
  c=colorbar
  set(c,'FontSize',FS)
  %    xlabel(c,'B_r (nT)')
  title('B_r (nT)')
  xlabel('UT')
  ylabel('magnetic lat. (deg.)')
  hold on;
  ax=axis;
  plot([t(1), t(end)],[mlatsrc, mlatsrc],'r--','MarkerSize',10,'LineWidth',2);
  hold off;
  
  subplot(132);
  param=Btheta*1e9;
  imagesc(t,mlat,param);
  colormap(parula(256));
  set(gca,'FontSize',FS);
  axis xy;
  datetick;
  axis tight;
  caxlim=max(abs(param(:)))
%  caxlim=max(caxlim,0.05)
  caxis([-caxlim,caxlim]);
  c=colorbar
  set(c,'FontSize',FS)
  %    xlabel(c,'B_\theta (nT)')
  title('B_\theta (nT)')
  xlabel('UT')
  ylabel('magnetic lat. (deg.)')
  hold on;
  ax=axis;
  plot([t(1), t(end)],[mlatsrc, mlatsrc],'r--','MarkerSize',10,'LineWidth',2);
  hold off;
  
  subplot(133);
  param=Bphi*1e9;
  imagesc(t,mlat,param);
  colormap(parula(256));
  set(gca,'FontSize',FS);
  axis xy;
  datetick;
  axis tight;
  caxlim=max(abs(param(:)))
%  caxlim=max(caxlim,0.05)
  caxis([-caxlim,caxlim]);
  c=colorbar
  set(c,'FontSize',FS)
  %    xlabel(c,'B_\phi (nT)')
  title('B_\phi (nT)')
  xlabel('UT')
  ylabel('magnetic lat. (deg.)')
  hold on;
  ax=axis;
  plot([t(1), t(end)],[mlatsrc, mlatsrc],'r--','MarkerSize',10,'LineWidth',2);
  hold off;
  
  
  print('-dpng',[direc,'/magplots/mag_timeseries.png'],'-r300');
end


%SAVE THE DATA TO A .MAT FILE IN CASE WE3 NEED IT LATER
t=datenum(simdate_series);
save([direc,'/magfields.mat'],'mlat','mlon','t','simdate_series','Brt','Bthetat','Bphit','-v7');


%RESET PATHS
rmpath ../script_utils;
