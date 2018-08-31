clear;
addpath ~/gemini/numerical;
simname='chile';
direc=['~/gemini/output/',simname];
direc2='~/gemini/output/chile_control';

load([direc,'/20150916_82473.mat']);
lx1=xg.lx(1); lx2=xg.lx(2);
lh=lx1;
lsp=size(ns,3);
xgsim=xg;


%NEW (PLOT) GRID IN R,TH
Re=6370e3;
%lth=750;
%lr=750;
lth=1000;
lr=1000;
rvals=xg.r(1:lh,:);
thvals=xg.theta(1:lh,:);
rmin=min(rvals(:));
rmax=max(rvals(:));
thmin=min(thvals(:));
thmax=max(thvals(:));
theta=linspace(thmin,thmax,lth);
r=linspace(rmin,rmax,lr)';
[THETA,R]=meshgrid(theta,r);
q=(Re./R).^2.*cos(THETA);
p=R./Re./sin(THETA).^2;


%OTHER GRID PARAMS FOR CALCULATING MAGNETIC FIELDS
zpr=R-6370e3;
xpr=THETA*(6370e3+100e3);
dzpr=zpr(2,1)-zpr(1,1);
dxpr=xpr(1,2)-xpr(1,1);

fprintf('(dz, dx)=(%d, %d)\n',dzpr,dxpr);

mlatrange=10;
izs=find(r>=100e3+6370e3 & r<=400e3+6370e3);
thmin=pi/2-pi/180*(29.212522+mlatrange);
thmax=pi/2-pi/180*(-29.212522-mlatrange);
ixs=find(theta>=thmin & theta<=thmax);
alt2=r(izs)/1e3-6370;
mlat2=90-fliplr(theta(ixs))*180/pi;

zmin=zpr(izs(1),1);
zmax=zpr(izs(end),1);
xmin=xpr(1,ixs(1));
xmax=xpr(1,ixs(end));

lint=500;
z=linspace(0,zmax,lint);
x=linspace(xmin,xmax,lint);
[x,z]=meshgrid(x,z);

mlatmin=mlat2(1);
mlatmax=mlat2(end);
mlat3=linspace(mlatmin,mlatmax,lint);
alt3=z/1e3;

lzpr=size(zpr,1);
lxpr=size(xpr,2);
lz=size(z,1);
lx=size(x,2);

mu0=4*pi*1e-7;


%DO INTERPOLATION IN Q,P SPACE SO WE CAN USE INTERP2
zvals=xg.z(1:lx1,:);
xvals=xg.x(1:lx1,:);
%itop=min(find(r>450e3+6370e3));
%times=20783:10:20783+3600;
times=82473:10:82473+3600;
lk=numel(times);
Bx=zeros(lk,lx);
By=zeros(lk,lx);
Bz=zeros(lk,lx);
simdateplot=[];
for k=1:lk
    %GET CURRENT DENSITIES FROM THIS TIME
    dmyplot=[16,9,2015];
    UTplot=(times(k)-1)/3600;
    loc=pwd;
    cd ~/gemini/vis;
    loadframe;
    cd(loc);

    simdateplotnow=[2015,09,16,times(k)/3600,0,0];
    simdateplot=[simdateplot;simdateplotnow];

    J1=Jx1(1:lh,:);
    J2=Jx2(1:lh,:);
    J3=Jx3(1:lh,:);


    %FIELD COORD INTERPS
    fprintf('Initiating interpolations for time frame %d...\n',it);
    J1I=interp2(xgsim.x2,xgsim.x1(1:lh),J1,p(:),q(:));
    J1I=reshape(J1I,size(R));
    inds=find(isnan(J1I));
    J1I(inds)=0;

    J2=Jx2(1:lh,:);
    J2I=interp2(xgsim.x2,xgsim.x1(1:lh),J2,p(:),q(:));
    J2I=reshape(J2I,size(R));
    inds=find(isnan(J2I));
    J2I(inds)=0;


    %J3 IS ALREADY GEOGRAPHIC ZONAL, BUT NEED TO ROTATE J1,2 AND INTEPOLATE
    Jr=J1.*dot(xg.er,xg.e1,3)+J2.*dot(xg.er,xg.e2,3);
    JrI=interp2(xgsim.x2,xgsim.x1(1:lh),Jr,p(:),q(:));
    JrI=reshape(JrI,size(R));
    inds=find(isnan(JrI));
    JrI(inds)=0;

    Jtheta=J1.*dot(xg.etheta,xg.e1,3)+J2.*dot(xg.etheta,xg.e2,3);
    JthetaI=interp2(xgsim.x2,xgsim.x1(1:lh),Jtheta,p(:),q(:));
    JthetaI=reshape(JthetaI,size(R));
    inds=find(isnan(JthetaI));
    JthetaI(inds)=0;

    Jphi=J3;
    JphiI=interp2(xgsim.x2,xgsim.x1(1:lh),Jphi,p(:),q(:));
    JphiI=reshape(JphiI,size(R));
    inds=find(isnan(JphiI));
    JphiI(inds)=0;


    %INTEGRATIONS
    fprintf('Beginning integrations for time step %d...\n',it);

    mu0=4*pi*1e-7;   %loadframe may be overwriting this with parallel mobility.  Yikes...


    %NOTE THAT Bx=Br (radial), By=Btheta (meridional) Bz=Bphi (zonal).  LIKEWISE FOR CURRENT DENSITY COMPONENTS.
    for iz=1:1
      for ix=1:lx
        Rz=z(iz,ix)-zpr;  %actually this is R1 (one-component)
        Rx=x(iz,ix)-xpr;  %R2
        dist2=max(Rz.^2+Rx.^2,dzpr.^2+dxpr.^2);

        integrand=mu0/4/pi*(-2*JphiI.*Rx)./dist2;
        intmids=zeros(lzpr-1,lxpr-1);
        intmids=1/4*integrand(1:lzpr-1,1:lxpr-1)+1/4*integrand(1:lzpr-1,2:lxpr)+1/4*integrand(2:lzpr,1:lxpr-1)+1/4*integrand(2:lzpr,2:lxpr);
        Bx(k,ix)=sum(intmids(:)*dzpr*dxpr);

        integrand=mu0/4/pi*(2*JphiI.*Rz)./dist2;
        intmids=zeros(lz-1,lx-1);
        intmids=1/4*integrand(1:lzpr-1,1:lxpr-1)+1/4*integrand(1:lzpr-1,2:lxpr)+1/4*integrand(2:lzpr,1:lxpr-1)+1/4*integrand(2:lzpr,2:lxpr); 
        By(k,ix)=sum(intmids(:)*dzpr*dxpr);

        integrand=mu0/4/pi*2*(JrI.*Rx-JthetaI.*Rz)./dist2;
        intmids=zeros(lz-1,lx-1);
        intmids=1/4*integrand(1:lzpr-1,1:lxpr-1)+1/4*integrand(1:lzpr-1,2:lxpr)+1/4*integrand(2:lzpr,1:lxpr-1)+1/4*integrand(2:lzpr,2:lxpr);
        Bz(k,ix)=sum(intmids(:)*dzpr*dxpr);
      end
    end

%{
    for iz=1:1   %only ground-level trace
      for ix=1:lx
        Rz=z(iz,ix)-zpr;
        Rx=x(iz,ix)-xpr;
        dist2=max(Rz.^2+Rx.^2,dzpr.^2+dxpr.^2);
  
        integrand=mu0/4/pi*2*JphiI.*Rz./dist2;
        intmids=zeros(lzpr-1,lxpr-1);
        intmids=1/4*integrand(1:lzpr-1,1:lxpr-1)+1/4*integrand(1:lzpr-1,2:lxpr)+1/4*integrand(2:lzpr,1:lxpr-1)+1/4*integrand(2:lzpr,2:lxpr);
        Bx(k,ix)=sum(intmids(:)*dzpr*dxpr);

        integrand=-1*mu0/4/pi*2*(JthetaI.*Rz-JrI.*Rx)./dist2;
        intmids=zeros(lz-1,lx-1);
        intmids=1/4*integrand(1:lzpr-1,1:lxpr-1)+1/4*integrand(1:lzpr-1,2:lxpr)+1/4*integrand(2:lzpr,1:lxpr-1)+1/4*integrand(2:lzpr,2:lxpr);
        By(k,ix)=sum(intmids(:)*dzpr*dxpr);

        integrand=mu0/4/pi*2*JphiI.*Rx./dist2;
        intmids=zeros(lz-1,lx-1);
        intmids=1/4*integrand(1:lzpr-1,1:lxpr-1)+1/4*integrand(1:lzpr-1,2:lxpr)+1/4*integrand(2:lzpr,1:lxpr-1)+1/4*integrand(2:lzpr,2:lxpr);
        Bz(k,ix)=sum(intmids(:)*dzpr*dxpr);
      end
    end
%}

end
save('-v7',['magfields_',simname,'_geogcomps_time.mat'],'Bx','By','Bz','mlat3','simdate');



%NOW DO THE PLOTS, FIELD COORDS. FIRST
fprintf('Initiating plotting...\n');
FS=18;

figure;
set(gcf,'PaperPosition',[0 0 8.5 8.5]);
tax=datenum(simdateplot);

subplot(311)
imagesc(tax,mlat3,1e9*fliplr(Bx)')
set(gca,'FontSize',FS);
datetick;
%datetick(15);
axis xy;
axis tight;
xlabel('UT')
ylabel('magnetic lat.')
title('vertical component of B')
% ix1s=find(r>=100e3+6370e3 & r<=400e3+6370e3);
Jrfilt=Bx;
maxJ=max(1e9*Jrfilt(:));
minJ=min(1e9*Jrfilt(:));
Jlim=max(abs(maxJ),abs(minJ));
caxis([-Jlim,Jlim]);
c=colorbar
set(c,'Fontsize',FS)
xlabel(c,'(nT)')
hold on;
plot(linspace(min(tax),max(tax),25),-20.41*ones(1,25),'w--','Linewidth',2);
axis tight;
hold off;


subplot(312)
imagesc(tax,mlat3,1e9*fliplr(By)')
set(gca,'FontSize',FS);
datetick;
%datetick(15);
axis xy;
axis tight;
xlabel('UT')
ylabel('magnetic lat.')
title('meridional component of B')
Jthfilt=By;
maxJ=max(1e9*Jthfilt(:));
minJ=min(1e9*Jthfilt(:));
Jlim=max(abs(maxJ),abs(minJ));
caxis([-Jlim,Jlim]);
c=colorbar
set(c,'Fontsize',FS)
xlabel(c,'(nT)')
hold on;
plot(linspace(min(tax),max(tax),25),-20.41*ones(1,25),'w--','Linewidth',2);
axis tight;
hold off;


subplot(313)
imagesc(tax,mlat3,1e9*fliplr(Bz)')
set(gca,'FontSize',FS);
datetick;
%datetick(15);
axis xy;
axis tight;
xlabel('UT')
ylabel('magnetic lat.')
title('zonal component of B')
Jphifilt=Bz;
maxJ=max(1e9*Jphifilt(:));
minJ=min(1e9*Jphifilt(:));
Jlim=max(abs(maxJ),abs(minJ));
Jlim=max(Jlim,0.1);
caxis([-Jlim,Jlim]);
c=colorbar
set(c,'Fontsize',FS)
xlabel(c,'(nT)')
hold on;
plot(linspace(min(tax),max(tax),25),-20.41*ones(1,25),'w--','Linewidth',2);
axis tight;
hold off;

print('-depsc2',['magfields_',simname,'_geogcomps_time.eps']);

