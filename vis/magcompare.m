%LOAD THE DATA
direc='./';
load([direc,'magfields.mat']);

[X,Y,Z]=ndgrid(x,y,z);
mu0=4*pi*1e-7;


%LOOP THROUGH A TAKE CURLS TO CHECK
[lx,ly,lz,lt]=size(Brt);
[lxp,lyp,lzp,lt]=size(Jrt);
Jxval=zeros(lx,ly,lz,lt);
Jyval=zeros(lx,ly,lz,lt);
Jzval=zeros(lx,ly,lz,lt);
for it=1:lt
  %STORE A TEMP. AUX. MAG FIELD
  Hx=Brt(:,:,:,it)/mu0;
  Hy=Bthetat(:,:,:,it)/mu0;
  Hz=Bphit(:,:,:,it)/mu0;
  %[Jxtmp,Jytmp,Jztmp]=curl(x,y,z,Hx,Hy,Hz);
  
  
  %TAKE THE CURL
  dHzdy=zeros(lx,ly,lz);
  dHzdy(:,1,:)=(Hz(:,2,:)-Hz(:,1,:))./(Y(:,2,:)-Y(:,1,:));
  dHzdy(:,2:end-1,:)=(Hz(:,3:end,:)-Hz(:,1:end-2,:))./(Y(:,3:end,:)-Y(:,1:end-2,:));
  dHzdy(:,end,:)=(Hz(:,end,:)-Hz(:,end-1,:))./(Y(:,end,:)-Y(:,end-1,:));
  
  dHydz=zeros(lx,ly,lz);
  dHydz(:,:,1)=(Hy(:,:,2)-Hy(:,:,1))./(Z(:,:,2)-Z(:,:,1));
  dHydz(:,:,2:end-1)=(Hy(:,:,3:end)-Hy(:,:,1:end-2))./(Z(:,:,3:end)-Z(:,:,1:end-2));
  dHydz(:,:,end)=(Hy(:,:,end)-Hy(:,:,end-1))./(Z(:,:,end)-Z(:,:,end-1));

  Jxtmp=dHzdy-dHydz;
  
  dHzdx=zeros(lx,ly,lz);
  dHzdx(1,:,:)=(Hz(2,:,:)-Hz(1,:,:))./(X(2,:,:)-X(1,:,:));
  dHzdx(2:end-1,:,:)=(Hz(3:end,:,:)-Hz(1:end-2,:,:))./(X(3:end,:,:)-X(1:end-2,:,:));
  dHzdx(end,:,:)=(Hz(end,:,:)-Hz(end-1,:,:))./(X(end,:,:)-X(end-1,:,:));
  
  dHxdz=zeros(lx,ly,lz);  
  dHxdz(:,:,1)=(Hx(:,:,2)-Hx(:,:,1))./(Z(:,:,2)-Z(:,:,1));
  dHxdz(:,:,2:end-1)=(Hx(:,:,3:end)-Hx(:,:,1:end-2))./(Z(:,:,3:end)-Z(:,:,1:end-2));
  dHxdz(:,:,end)=(Hx(:,:,end)-Hx(:,:,end-1))./(Z(:,:,end)-Z(:,:,end-1));
  
  Jytmp=-1*(dHzdx-dHxdz);
  
  dHydx=zeros(lx,ly,lz);
  dHydx(1,:,:)=(Hy(2,:,:)-Hy(1,:,:))./(X(2,:,:)-X(1,:,:));
  dHydx(2:end-1,:,:)=(Hy(3:end,:,:)-Hy(1:end-2,:,:))./(X(3:end,:,:)-X(1:end-2,:,:));
  dHydx(end,:,:)=(Hy(end,:,:)-Hy(end-1,:,:))./(X(end,:,:)-X(end-1,:,:));
  
  dHxdy=zeros(lx,ly,lz);
  dHxdy(:,1,:)=(Hx(:,2,:)-Hx(:,1,:))./(Y(:,2,:)-Y(:,1,:));
  dHxdy(:,2:end-1,:)=(Hx(:,3:end,:)-Hx(:,1:end-2,:))./(Y(:,3:end,:)-Y(:,1:end-2,:));
  dHxdy(:,end,:)=(Hx(:,end,:)-Hx(:,end-1,:))./(Y(:,end,:)-Y(:,end-1,:));  
  
  Jztmp=dHydx-dHxdy;
  
  
  %STORE THE RESULTS
  Jxval(:,:,:,it)=Jxtmp;
  Jyval(:,:,:,it)=Jytmp;
  Jzval(:,:,:,it)=Jztmp;  
end


%SOME PLOTS TO CHECK
figure;
title('x-y slice')

iphi=floor(lz/2);
zref=z(iphi);
[tmp,iphip]=min(abs(zp-zref));
%iphip=floor(lzp/2);
it=1;

subplot(231);
imagesc(y,x,Jxval(:,:,iphi,it));
axis xy;
ax=axis;
colorbar;

subplot(232);
imagesc(y,x,Jyval(:,:,iphi,it));
axis xy;
colorbar;

subplot(233);
imagesc(y,x,Jzval(:,:,iphi,it));
axis xy;
colorbar;

subplot(234);
imagesc(yp,xp,Jrt(:,:,iphip,it));
axis xy;
axis(ax);
colorbar;

subplot(235);
imagesc(yp,xp,Jthetat(:,:,iphip,it));
axis xy;
axis(ax);
colorbar;

subplot(236);
imagesc(yp,xp,Jphit(:,:,iphip,it));
axis xy;
axis(ax);
colorbar;


figure;
title('x-z slice')

itheta=floor(ly/2);
yref=y(itheta);
[tmp,ithetap]=min(abs(yp-yref));
it=1;

subplot(231);
imagesc(z,x,squeeze(Jxval(:,itheta,:,it)));
axis xy;
ax=axis;
colorbar;

subplot(232);
imagesc(z,x,squeeze(Jyval(:,itheta,:,it)));
axis xy;
colorbar;

subplot(233);
imagesc(z,x,squeeze(Jzval(:,itheta,:,it)));
axis xy;
colorbar;

subplot(234);
imagesc(zp,xp,squeeze(Jrt(:,ithetap,:,it)));
axis xy;
axis(ax);
colorbar;

subplot(235);
imagesc(zp,xp,squeeze(Jthetat(:,ithetap,:,it)));
axis xy;
axis(ax);
colorbar;

subplot(236);
imagesc(zp,xp,squeeze(Jphit(:,ithetap,:,it)));
axis xy;
axis(ax);
colorbar;


figure;
title('y-z slice')

ir=floor(lx/2);
xref=x(ir);
[tmp,irp]=min(abs(xp-xref));
it=1;

subplot(231);
imagesc(y,z,squeeze(Jxval(ir,:,:,it))');
axis xy;
ax=axis;
colorbar;

subplot(232);
imagesc(y,z,squeeze(Jyval(ir,:,:,it))');
axis xy;
colorbar;

subplot(233);
imagesc(y,z,squeeze(Jzval(ir,:,:,it))');
axis xy;
colorbar;

subplot(234);
imagesc(yp,zp,squeeze(Jrt(irp,:,:,it))');
axis xy;
axis(ax);
colorbar;

subplot(235);
imagesc(yp,zp,squeeze(Jthetat(irp,:,:,it))');
axis xy;
axis(ax);
colorbar;

subplot(236);
imagesc(yp,zp,squeeze(Jphit(irp,:,:,it))');
axis xy;
axis(ax);
colorbar;