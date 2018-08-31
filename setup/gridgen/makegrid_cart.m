function xg=makegrid_cart(xmin,xmax,I,lperp,glat,glon)


%SETUP NONUNIFORM GRID IN ALTITUDE AND FIELD LINE DISTANCE
altmin=80e3;
altmax=1000e3;
%alt=linspace(altmin,altmax,920);
ialt=1;
alt(ialt)=altmin;
while alt(ialt)<altmax
    ialt=ialt+1;
    dalt=10e3+8e3*tanh((alt(ialt-1)-500e3)/150e3);
    alt(ialt)=alt(ialt-1)+dalt;
end
% while alt(ialt)<altmax
%     ialt=ialt+1;
%     dalt=10+9.5*tanh((alt(ialt-1)-500)/150);
%     alt(ialt)=alt(ialt-1)+dalt;
% end
alt=alt(:);
z=alt*cscd(I);
lz=numel(z);
lx1=lz;


%TRANSVERSE GRID (BASED ON SIZE OF CURRENT REGION SPECIFIED ABOVE)
x=linspace(xmin,xmax,lperp);
lx=lperp;
lx2=lx;


%COMPUTE CELL WALL LOCATIONS
xi=zeros(1,lx+1);
xi(2:lx)=1/2*(x(2:lx)+x(1:lx-1));
xi(1)=x(1)-1/2*(x(2)-x(1));
xi(lx+1)=x(lx)+1/2*(x(lx)-x(lx-1));
zi=zeros(lz+1,1);
zi(2:lz)=1/2*(z(2:lz)+z(1:lz-1));
zi(1)=z(1)-1/2*(z(2)-z(1));
zi(lz+1)=z(lz)+1/2*(z(lz)-z(lz-1));


%GRAVITATIONAL FIELD COMPONENTS IN DIPOLE SYSTEM
Re=6370e3;
G=6.67428e-11;
Me=5.9722e24;
r=alt+Re;
g=G*Me./r.^2;
gz=repmat(-1*g,1,lperp);


%STORE RESULTS IN GRID DATA STRUCTURE
xg.x1=z; xg.x2=x; xg.x3=0; 
xg.x1i=zi; xg.x2i=xi;
lx=[numel(xg.x1),numel(xg.x2)];
xg.lx=lx;

dx1=xg.x1(2)-xg.x1(1);
dxn=xg.x1(lx(1))-xg.x1(lx(1)-1);
xg.dx1f=[xg.x1(2:lx(1))-xg.x1(1:lx(1)-1); dxn];         %FWD DIFF
xg.dx1b=[dx1; xg.x1(2:lx(1))-xg.x1(1:lx(1)-1)];         %BACK DIFF
xg.dx1h=xg.x1i(2:lx(1)+1)-xg.x1i(1:lx(1));              %MIDPOINT DIFFS

dx1=xg.x2(2)-xg.x2(1);
dxn=xg.x2(lx(2))-xg.x2(lx(2)-1);
xg.dx2f=[xg.x2(2:lx(2))-xg.x2(1:lx(2)-1), dxn];         %FWD DIFF
xg.dx2b=[dx1, xg.x2(2:lx(2))-xg.x2(1:lx(2)-1)];         %BACK DIFF
xg.dx2h=xg.x2i(2:lx(2)+1)-xg.x2i(1:lx(2));              %MIDPOINT DIFFS

xg.h1=ones(xg.lx); xg.h2=ones(xg.lx); xg.h3=ones(xg.lx);
xg.h1x1i=ones(lx(1)+1,lx(2)); xg.h2x1i=ones(lx(1)+1,lx(2)); xg.h3x1i=ones(lx(1)+1,lx(2));
xg.h1x2i=ones(lx(1),lx(2)+1); xg.h2x2i=ones(lx(1),lx(2)+1); xg.h3x2i=ones(lx(1),lx(2)+1);

xg.e1=[]; xg.e2=[]; xg.e3=[];

xg.r=repmat(r,1,lx(2)); xg.theta=[]; xg.phi=[];
xg.rx1i=[]; xg.thetax1i=[];
xg.rx2i=[]; xg.thetax2i=[];

xg.er=[]; xg.etheta=[]; xg.ephi=[];

xg.I=I*ones(1,lx2);

xg.x=x; xg.z=z;
xg.alt=xg.r-Re;

xg.gx1=gz; xg.gx2=zeros(xg.lx);

xg.Bmag=50000e-9;

xg.glat=glat*ones(lx1,lx2); xg.glon=glon*ones(lx1,lx2);

xg.xp=x; xg.zp=z;

xg.inull=[];
%xg.nullpts=[];
xg.nullpts=zeros(xg.lx(1),xg.lx(2));

end
