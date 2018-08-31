function h=plotslice3D_curv_corner_HTMS(t,xg,parm)

%CLEAR AND SET FIGURE HANDLES
%clf;
h=gcf;


%SIZE OF SIMULATION
lx1=xg.lx(1); lx2=xg.lx(2); lx3=xg.lx(3);
%inds1=3:lx1+2;
%inds2=3:lx2+2;
%inds3=3:lx3+2;
inds1=1:lx1;
inds2=1:lx2;
inds3=1:lx3;
lsp=size(parm,3);
Re=6370e3;


%SIZE OF PLOT GRID THAT WE ARE INTERPOLATING ONTO
izs=1:lx1;
iys=1:lx3;
ixs=1:lx2;
meantheta=mean(xg.theta(:));
meanphi=mean(xg.phi(:));
x=(xg.theta-meantheta);   %this is a mag colat. coordinate and is only used for defining grid in linspaces below
y=(xg.phi-meanphi);       %mag. lon coordinate
z=xg.alt/1e3;
lxp=500;
lyp=250;
lzp=500;
minx=min(x(:));
maxx=max(x(:));
miny=min(y(:));
maxy=max(y(:));
minz=min(z(:));
maxz=max(z(:));
xp=linspace(minx,maxx,lxp);
yp=linspace(miny,maxy,lyp);
zp=linspace(minz,maxz,lzp)';


%INTERPOLATE ONTO PLOTTING GRID
[X,Y,Z]=meshgrid(xp,yp,zp*1e3);

%{
    %WE NEED THE Q,P PAIRS FOR THIS GRID
    %BACK TO SPHERICAL FIRST (ROTATE)
    r1=mean(xg.r(1,:));
    r2=mean(xg.r(lx1,:));
    if abs(r1-r2)>10e3
        if r1<r2
            meanth=mean(xg.theta(1,:));
        else
            meanth=mean(xg.theta(lx1,:));
        end
    else
        meanth=pi/2;
    end
    rotmat=[cos(meanth), -sin(meanth); sin(meanth), cos(meanth)]';
    xx=X(:);
    zx=Z(:);
    xxp=zeros(size(xx));
    zxp=xxp;
    lk=numel(xx);
    for k=1:lk
        xz=rotmat*[xx(k); zx(k)];
        xxp(k)=xz(1); zxp(k)=xz(2);
    end
    
    %NOW OFFSET
    xsurf=Re*sin(meanth);
    zsurf=Re*cos(meanth);
    xxp=xxp+xsurf;
    zxp=zxp+zsurf;
    
    
    %NOW TO SPHERICAL
    rxp=sqrt(xxp.^2+zxp.^2);
    thetaxp=acos(zxp./rxp);
%}

%DIRECT TO SPHERICAL
rxp=Z(:)+Re;
thetaxp=X(:)+meantheta;
phixp=Y(:)+meanphi;


%NOW SPHERICAL TO DIPOLE
qplot=(Re./rxp).^2.*cos(thetaxp);
pplot=rxp/Re./sin(thetaxp).^2;
phiplot=phixp;    %phi is same in spherical and dipole


%NOW WE CAN DO A `PLAID' INTERPOLATION
parmtmp=permute(parm,[3,2,1]);
[X2,X3,X1]=meshgrid(xg.x2(inds2),xg.x3(inds3),xg.x1(inds1));
parmp=interp3(X2,X3,X1,parmtmp,pplot,phiplot,qplot);
parmp=reshape(parmp,lyp,lxp,lzp);    %slice expects the first dim. to be "y"


%CONVERT ANGULAR COORDINATES TO MLAT,MLON
xp=90-(xp+mean(xg.theta(:)))*180/pi;
[xp,inds]=sort(xp);
parmp=parmp(:,inds,:);

yp=(yp+mean(xg.phi(:)))*180/pi;
[yp,inds]=sort(yp);
parmp=parmp(inds,:,:);


%NOW THAT WE'VE SORTED, WE NEED TO REGENERATE THE MESHGRID
[XP,YP,ZP]=meshgrid(xp,yp,zp);


%MAKE THE PLOTS!
FS=16;
% %xslice=[90-meantheta*180/pi];   %make a few slices at differen altitudes
% xslice=[];
% yslice=[];
% zslice=[150,300,750];
% h=slice(XP,YP,ZP,parmp,xslice,yslice,zslice);
% shading flat;
hold on;
%ix3=floor(lx3/2);
ix3=lx3-4;    %don't pick the last one (in case it is not there)
surfx=90-(x(:,:,ix3)+meantheta)*180/pi;    %create a surface out of a meridional slice
surfy=(y(:,:,ix3)+meanphi)*180/pi;
surfz=z(:,:,ix3);
h=slice(XP,YP,ZP,parmp,surfx,surfy,surfz);
shading interp;
% %ix2=floor(lx2/2);
% ix2=2;
% surfx=squeeze(90-(x(:,ix2,:)+meantheta)*180/pi);    %create a another surface out of slice along constant L-shell
% surfy=squeeze((y(:,ix2,:)+meanphi)*180/pi);
% surfz=squeeze(z(:,ix2,:));
% h=slice(XP,YP,ZP,parmp,surfx,surfy,surfz);
% shading interp;
hold off;
%set(h,'alphadata',~isnan(parmp));
set(gca,'FontSize',FS);
axis xy;
axis tight;
%c=colorbar;
%xlabel(c,sprintf('n_{e} (log_{10} m^{-3})'))
%xlabel(c,sprintf('v_{i,||} (m/s)'))
%xlabel('mag. lat. (deg.)');
%ylabel('mag. long. (deg.)');
%zlabel('alt. (km)')


%CONSTRUCT A STRING FOR THE TIME AND DATE
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
  UTsecstr=['0',UTsecstr]
end

dmy=[0,0,0];
timestr=[UThrsstr,':',UTminstr,':',UTsecstr];
%strval=sprintf('%s \n %s',[num2str(dmy(1)),'/',num2str(dmy(2)),'/',num2str(dmy(3))], ...
%    [num2str(t),' UT']);
strval=sprintf('%s \n %s',[num2str(dmy(1)),'/',num2str(dmy(2)),'/',num2str(dmy(3))], ...
    [timestr,' UT']);
%text(xp(round(lxp/10)),zp(lzp-round(lzp/7.5)),strval,'FontSize',18,'Color',[0.66 0.66 0.66],'FontWeight','bold');
text(xp(round(lxp/10)),zp(lzp-round(lzp/7.5)),strval,'FontSize',16,'Color',[1 1 1],'FontWeight','bold');

end
