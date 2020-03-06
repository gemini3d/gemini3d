function h=plot3D_curv_frames_long(ymd,UTsec,xg,parm,parmlbl,caxlims,sourceloc,hf, cmap)

narginchk(4,9)
validateattr(ymd, {'numeric'}, {'vector', 'numel', 3}, mfilename, 'year month day', 1)
validateattr(UTsec, {'numeric'}, {'scalar'}, mfilename, 'UTC second', 2)
validateattr(xg, {'struct'}, {'scalar'}, mfilename, 'grid structure', 3)
validateattr(parm, {'numeric'}, {'real'}, mfilename, 'parameter to plot',4)
if nargin<5, parmlbl=''; end
validateattr(parmlbl, {'char'}, {'vector'}, mfilename, 'parameter label', 5)

if nargin<6
  caxlims=[];
else
  validateattr(caxlims, {'numeric'}, {'vector', 'numel', 2}, mfilename, 'plot intensity (min, max)', 6)
end

if nargin<7  || isempty(sourceloc) % leave || for validate
  sourceloc = [];
else
  validateattr(sourceloc, {'numeric'}, {'vector', 'numel', 2}, mfilename, 'source magnetic coordinates', 7)
end
if nargin<8 || isempty(hf)
  hf = figure();
end
if nargin<9 || isempty(cmap)
  cmap = parula(256);
end


%SOURCE LOCATION (SHOULD PROBABLY BE AN INPUT)
if (~isempty(sourceloc))
  sourcemlat=sourceloc(1);
  sourcemlon=sourceloc(2);
  flagsource=1;
else
  flagsource=0;
end

%SIZE OF SIMULATION
lx1=xg.lx(1); lx2=xg.lx(2); lx3=xg.lx(3);
inds1=3:lx1+2;
inds2=3:lx2+2;
inds3=3:lx3+2;
Re=6370e3;


%JUST PICK AN X3 LOCATION FOR THE MERIDIONAL SLICE PLOT, AND AN ALTITUDE FOR THE LAT./LON. SLICE
ix3=floor(lx3/2);
altref=300;


%SIZE OF PLOT GRID THAT WE ARE INTERPOLATING ONTO
meantheta=mean(xg.theta(:));
meanphi=mean(xg.phi(:));
x=(xg.theta-meantheta);   %this is a mag colat. coordinate and is only used for defining grid in linspaces below
y=(xg.phi-meanphi);       %mag. lon coordinate
z=xg.alt/1e3;
%lxp=500;
lxp=1500;
lyp=500;
lzp=500;
minx=min(x(:));
maxx=max(x(:));
miny=min(y(:));
maxy=max(y(:));
%minz=min(z(:));
minz=-10e0;     %to give some space for the marker on th plots
maxz=max(z(:));
xp=linspace(minx,maxx,lxp);
yp=linspace(miny,maxy,lyp);
zp=linspace(minz,maxz,lzp)';

%{
%ix1s=floor(lx1/2):lx1;    %only valide for a grid which is symmetric aboutu magnetic equator... (I think)
ix1s=find(xg.x1(inds1)>=0);    %works for asymmetric grids
minz=0;
maxz=max(xg.alt(:));
[tmp,ix1]=min(abs(xg.alt(ix1s,1,1)-maxz*1e3));
ix1=ix1s(ix1);
thetavals=xg.theta(ix1:lx1,:,:);
meantheta=mean(thetavals(:));
phivals=xg.phi(ix1:lx1,:,:);
meanphi=mean(phivals(:));
x=(thetavals-meantheta);      %this is a mag colat. coordinate and is only used for defining grid in linspaces below and the parametric surfaces in the plots
y=(phivals-meanphi);          %mag. lon coordinate
z=xg.alt(ix1:lx1,:,:)/1e3;    %altitude
lxp=500;
lyp=500;
lzp=500;
minx=min(x(:));
maxx=max(x(:));%+0.5*(max(x(:))-min(x(:)));
miny=min(y(:));
maxy=max(y(:));
xp=linspace(minx,maxx,lxp);
yp=linspace(miny,maxy,lyp);
zp=linspace(minz,maxz,lzp)';
%}


%INTERPOLATE ONTO PLOTTING GRID
[X,Z]=meshgrid(xp,zp*1e3);    %meridional meshgrid


%DIRECT TO SPHERICAL
rxp=Z(:)+Re;
thetaxp=X(:)+meantheta;
%phixp=Y(:)+meanphi;


%NOW SPHERICAL TO DIPOLE
qplot=(Re./rxp).^2.*cos(thetaxp);
pplot=rxp/Re./sin(thetaxp).^2;
%phiplot=phixp;    %phi is same in spherical and dipole


%-----NOW WE CAN DO A `PLAID' INTERPOLATION - THIS ONE IS FOR THE MERIDIONAL SLICE
parmtmp=parm(:,:,ix3);
parmp=interp2(xg.x2(inds2),xg.x1(inds1),parmtmp,pplot,qplot);
parmp=reshape(parmp,lzp,lxp);    %slice expects the first dim. to be "y" ("z" in the 2D case)


%LAT./LONG. SLICE COORDIANTES
zp2=[290,300,310];
lzp2=3;
[X2,Y2,Z2]=meshgrid(xp,yp,zp2*1e3);       %lat./lon. meshgrid, need 3D since and altitude slice cuts through all 3 dipole dimensions

rxp2=Z2(:)+Re;
thetaxp2=X2(:)+meantheta;
phixp2=Y2(:)+meanphi;

qplot2=(Re./rxp2).^2.*cos(thetaxp2);
pplot2=rxp2/Re./sin(thetaxp2).^2;
phiplot2=phixp2;    %phi is same in spherical and dipole


%-----NOW WE CAN DO A `PLAID' INTERPOLATION - THIS ONE IS FOR THE LAT/LON SLICE
parmtmp=permute(parm,[3,2,1]);
x3interp=xg.x3(inds3);
x3interp=x3interp(:);     %interp doesn't like it unless this is a column vector
parmp2=interp3(xg.x2(inds2),x3interp,xg.x1(inds1),parmtmp,pplot2,phiplot2,qplot2);
parmp2=reshape(parmp2,lyp,lxp,lzp2);    %slice expects the first dim. to be "y"


%ALT./LONG. SLICE COORDIANTES
if (flagsource)
  sourcetheta=pi/2-sourcemlat*pi/180;
else
  thetaref=mean(mean(xg.theta(1:floor(end/2),:,:)));
  sourcetheta=thetaref;
end
sourcex=sourcetheta-meantheta;
xp3=[sourcex-1*pi/180,sourcex,sourcex+1*pi/180];    %place in the plane of the source
lxp3=numel(xp3);
zp3=linspace(minz,750,500);
%minzp3=min(zp3);
%maxzp3=max(zp3);
lzp3=numel(zp3);
[X3,Y3,Z3]=meshgrid(xp3,yp,zp3*1e3);       %lat./lon. meshgrid, need 3D since and altitude slice cuts through all 3 dipole dimensions

rxp3=Z3(:)+Re;
thetaxp3=X3(:)+meantheta;
phixp3=Y3(:)+meanphi;

qplot3=(Re./rxp3).^2.*cos(thetaxp3);
pplot3=rxp3/Re./sin(thetaxp3).^2;
phiplot3=phixp3;    %phi is same in spherical and dipole


%-----NOW WE CAN DO A `PLAID' INTERPOLATION - THIS ONE IS FOR THE LONG/ALT SLICE
%parmtmp=permute(parm,[3,2,1]);
parmp3=interp3(xg.x2(inds2),x3interp,xg.x1(inds1),parmtmp,pplot3,phiplot3,qplot3);
parmp3=reshape(parmp3,lyp,lxp3,lzp3);    %slice expects the first dim. to be "y"


%CONVERT ANGULAR COORDINATES TO MLAT,MLON
xp=90-(xp+meantheta)*180/pi;
[xp,inds]=sort(xp);
parmp=parmp(:,inds);
parmp2=parmp2(:,inds,:);
%parmp3=parmp3(:,inds,:);

yp=(yp+meanphi)*180/pi;
[yp,inds]=sort(yp);
%parmp=parmp(inds,:,:);
parmp2=parmp2(inds,:,:);
parmp3=parmp3(inds,:,:);


%COMPUTE SOME BOUNDS FOR THE PLOTTING
minxp=min(xp(:));
maxxp=max(xp(:));
%minyp=min(yp(:));
%maxyp=max(yp(:));
minzp=min(zp(:));
maxzp=max(zp(:));


%NOW THAT WE'VE SORTED, WE NEED TO REGENERATE THE MESHGRID
%[XP,YP,ZP]=meshgrid(xp,yp,zp);
FS=8;

%MAKE THE PLOT!
ha=subplot(1,3,1, 'parent', hf, 'nextplot', 'add', 'FontSize',FS);
h=imagesc(ha,xp,zp,parmp);
if (flagsource)
  plot(ha,[minxp,maxxp],[altref,altref],'--','LineWidth',1,'Color',[0.25 0.25 0.25]);
  plot(ha,[sourcemlat,sourcemlat],[minzp,maxzp],'k--','LineWidth',1);
  plot(ha,sourcemlat,0,'r^','MarkerSize',6,'LineWidth',2);
end
set(h,'alphadata',~isnan(parmp));

tight_axis(ha)
colormap(ha,cmap)
caxis(ha,caxlims);
c=colorbar(ha);
xlabel(c,parmlbl);
xlabel(ha,'magnetic latitude (deg.)');
ylabel(ha,'altitude (km)');



%CONSTRUCT A STRING FOR THE TIME AND DATE
t=UTsec/3600;     %decimal hours
UThrs=floor(t);
UTmin=floor((t-UThrs)*60);
UTsec=floor((t-UThrs-UTmin/60)*3600);
UThrsstr=num2str(UThrs);
UTminstr=num2str(UTmin);
if (numel(UTminstr)==1)
  UTminstr=['0',UTminstr];
end
UTsecstr=num2str(UTsec);
if (numel(UTsecstr)==1)
  UTsecstr=['0',UTsecstr];
end

timestr=[UThrsstr,':',UTminstr,':',UTsecstr];
strval=sprintf('%s \n %s',[num2str(ymd(1)),'/',num2str(ymd(2)),'/',num2str(ymd(3))], ...
    [timestr,' UT']);
title(ha,strval,'Color','black');



ha=subplot(1,3,2, 'parent', hf, 'nextplot', 'add', 'FontSize',FS);
%h=imagesc(xp,yp,parmp2(:,:,2));
h=imagesc(ha,yp,xp,parmp2(:,:,2)');    %so latitude is the vertical axis
if (flagsource)
  %plot([minxp,maxxp],[sourcemlon,sourcemlon],'w--','LineWidth',2);
  %plot(sourcemlat,sourcemlon,'r^','MarkerSize',12,'LineWidth',2);
%  plot([sourcemlon,sourcemlon],[minxp,maxxp],'w--','LineWidth',1);
  plot(ha,sourcemlon,sourcemlat,'r^','MarkerSize',6,'LineWidth',2);
end
%set(h,'alphadata',~isnan(parmp2(:,:,2)));
set(h,'alphadata',~isnan(parmp2(:,:,2)'));

tight_axis(ha)
colormap(ha,cmap)
caxis(ha,caxlims);
c=colorbar(ha);
xlabel(c,parmlbl);
%xlabel('magnetic latitude (deg.)');
%ylabel('magnetic longitude (deg.)');
xlabel(ha,'magnetic long. (deg.)');
ylabel(ha,'magnetic lat. (deg.)');


ha=subplot(1,3,3, 'parent', hf, 'nextplot', 'add', 'FontSize',FS);
h=imagesc(ha,yp,zp3,squeeze(parmp3(:,2,:))');    %so latitude is the vertical axis
if (flagsource)
%  plot([sourcemlon,sourcemlon],[minzp3,maxzp3],'k--','LineWidth',1);
  plot(ha,sourcemlon,0,'r^','MarkerSize',6,'LineWidth',2);
end
%set(h,'alphadata',~isnan(squeeze(parmp3(:,2,:))'));

tight_axis(ha)
colormap(ha,cmap)
caxis(ha,caxlims);
c=colorbar(ha);
xlabel(c,parmlbl);
xlabel(ha,'magnetic long. (deg.)');
ylabel(ha,'altitude (km)');


%CONSTRUCT A STRING FOR THE TIME AND DATE
%ha=subplot(1,3,1);

%ttxt = time2str(ymd, UTsec);
%title(ha, ttxt)
%%text(xp(round(lxp/10)),zp(lzp-round(lzp/7.5)),strval,'FontSize',18,'Color',[0.66 0.66 0.66],'FontWeight','bold');
%%text(xp(round(lxp/10)),zp(lzp-round(lzp/7.5)),strval,'FontSize',16,'Color',[0.5 0.5 0.5],'FontWeight','bold');

end
