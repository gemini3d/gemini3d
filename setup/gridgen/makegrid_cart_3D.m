function xgf = makegrid_cart_3D(p)

narginchk(1, 1)
validateattributes(p, {'struct'}, {'scalar'})

%% ADD IN GHOST CELLS (Z WILL BE HANDLED LATER)
lx = p.lxp + 4;
ly = p.lyp + 4;

%% create altitude grid
%p.alt_min = 80e3;
%p.alt_max = 1000e3;
%p.alt_scale = [10e3, 8e3, 500e3, 150e3];
if isfield(p, 'alt_min') && isfield(p, 'alt_max') && isfield(p, 'alt_scale') && isfield(p,'I')
  z = altitude_grid(p.alt_min, p.alt_max, p.I, p.alt_scale);
else
  disp(['makegrid_cart_3D: reusing grid from ', p.eqdir])
  xeq = readgrid(p.eqdir, p.format);
  z = xeq.x1;
  clear('xeq')
end
%% TRANSVERSE GRID (BASED ON SIZE OF CURRENT REGION SPECIFIED ABOVE)
%EAST
xmin = -p.xdist/2;
xmax = p.xdist/2;

if p.lxp == 1  % degenerate dimension
  x = linspace(xmin, xmax, lx);
else
  % exclude the ghost cells when setting extents
  x = linspace(xmin, xmax, p.lxp);
  dx1 = x(2) - x(1);
  dxn = x(end) - x(end-1);
  % now tack on ghost cells so they are outside user-specified region
  x=[x(1)-2*dx1, x(1)-dx1,x, x(end)+dxn, x(end)+2*dxn];
end
lx2=lx;

% NORTH
ymin = -p.ydist/2;
ymax = p.ydist/2;

if p.lyp == 1  % degenerate dimension
  y = linspace(ymin, ymax, ly);
else
  % exclude the ghost cells when setting extents
  y = linspace(ymin, ymax, p.lyp);
  dy1=y(2)-y(1);
  dyn=y(end)-y(end-1);
  % now tack on ghost cells so they are outside user-specified region
  y=[y(1)-2*dy1,y(1)-dy1,y,y(end)+dyn,y(end)+2*dyn];
end
lx3=ly;

%% COMPUTE CELL WALL LOCATIONS
xi=zeros(1,lx+1);
xi(2:lx)=1/2*(x(2:lx) + x(1:lx-1));
xi(1)=x(1)-1/2 * (x(2) - x(1));
xi(lx+1)=x(lx)+1/2*(x(lx)-x(lx-1));

yi=zeros(1,ly+1);
yi(2:ly)=1/2*(y(2:ly)+y(1:ly-1));
yi(1)=y(1)-1/2*(y(2)-y(1));
yi(ly+1)=y(ly)+1/2*(y(ly)-y(ly-1));

lx1 = length(z);
zi=zeros(lx1+1,1);
zi(2:lx1)=1/2*(z(2:lx1)+z(1:lx1-1));
zi(1)=z(1)-1/2*(z(2)-z(1));
zi(lx1+1)=z(lx1)+1/2*(z(lx1)-z(lx1-1));


%% GRAVITATIONAL FIELD COMPONENTS IN DIPOLE SYSTEM
Re=6370e3;
G=6.67428e-11;
Me=5.9722e24;
r=z+Re;
g=G*Me./r.^2;
gz=repmat(-1*g,[1,lx2,lx3]);


%DISTANCE EW AND NS (FROM ENU (or UEN in our case - cyclic permuted) COORD. SYSTEM) NEED TO BE CONVERTED TO DIPOLE SPHERICAL AND THEN
%GLAT/GLONG - BASICALLY HERE WE ARE MAPPING THE CARTESIAN GRID ONTO THE
%SURFACE OF A SPHERE THEN CONVERTING TO GEOGRAPHIC.
[thetactr,phictr] = geog2geomag(p.glat, p.glon);    %get the magnetic coordinates of the grid center, based on user input

%% Center of earth distance
r=Re+z;
r=reshape(r,[lx1,1,1]);
r=repmat(r(:),[1,lx2,lx3]);

%% Northward angular distance
gamma2=y/Re;    %must retain the sign of x3
theta=thetactr-gamma2;   %minus because distance north is against theta's direction
theta=reshape(theta,[1,1,lx3]);
theta=repmat(theta,[lx1,lx2,1]);

%% Eastward angular distance
%gamma1=x/Re;     %must retain the sign of x2
gamma1=x/Re/sin(thetactr);     %must retain the sign of x2, just use theta of center of grid
phi=phictr+gamma1;
phi=reshape(phi,[1,lx2,1]);
phi=repmat(phi,[lx1,1,lx3]);

%% COMPUTE THE GEOGRAPHIC COORDINATES OF EACH GRID POINT
[glatgrid,glongrid]=geomag2geog(theta,phi);

%% COMPUTE ECEF CARTESIAN IN CASE THEY ARE NEEDED
xECEF=r.*sin(theta).*cos(phi);
yECEF=r.*sin(theta).*sin(phi);
zECEF=r.*cos(theta);

%% COMPUTE SPHERICAL ECEF UNIT VECTORS - CARTESIAN-ECEF COMPONENTS
disp('MAKEGRID_CART_3D.M --> Calculating spherical unit vectors')
er(:,:,:,1)=sin(theta).*cos(phi);    %xECEF-component of er
er(:,:,:,2)=sin(theta).*sin(phi);    %yECEF
er(:,:,:,3)=cos(theta);              %zECEF
etheta(:,:,:,1)=cos(theta).*cos(phi);
etheta(:,:,:,2)=cos(theta).*sin(phi);
etheta(:,:,:,3)=-sin(theta);
ephi(:,:,:,1)=-sin(phi);
ephi(:,:,:,2)=cos(phi);
ephi(:,:,:,3)=zeros(lx1,lx2,lx3);


%UEN UNIT VECTORS IN ECEF COMPONENTS
e1=er;    %up is the same cirection as from ctr of earth
e2=ephi;    %e2 is same as ephi
e3=-1*etheta;    %etheta is positive south, e3 is pos. north


%STORE RESULTS IN GRID DATA STRUCTURE
xg.x1=z; xg.x2=x; xg.x3=y;
xg.x1i=zi; xg.x2i=xi; xg.x3i=yi;
lx=[numel(xg.x1),numel(xg.x2),numel(xg.x3)];
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

dx1=xg.x3(2)-xg.x3(1);
dxn=xg.x3(lx(3))-xg.x3(lx(3)-1);
xg.dx3f=[xg.x3(2:lx(3))-xg.x3(1:lx(3)-1), dxn];         %FWD DIFF
xg.dx3b=[dx1, xg.x3(2:lx(3))-xg.x3(1:lx(3)-1)];         %BACK DIFF
xg.dx3h=xg.x3i(2:lx(3)+1)-xg.x3i(1:lx(3));              %MIDPOINT DIFFS

xg.h1=ones(xg.lx); xg.h2=ones(xg.lx); xg.h3=ones(xg.lx);
xg.h1x1i=ones(lx(1)+1,lx(2),lx(3)); xg.h2x1i=ones(lx(1)+1,lx(2),lx(3)); xg.h3x1i=ones(lx(1)+1,lx(2),lx(3));
xg.h1x2i=ones(lx(1),lx(2)+1,lx(3)); xg.h2x2i=ones(lx(1),lx(2)+1,lx(3)); xg.h3x2i=ones(lx(1),lx(2)+1,lx(3));
xg.h1x3i=ones(lx(1),lx(2),lx(3)+1); xg.h2x3i=ones(lx(1),lx(2),lx(3)+1); xg.h3x3i=ones(lx(1),lx(2),lx(3)+1);

%Cartesian, ECEF representation of curvilinar coordinates
xg.e1=e1; xg.e2=e2; xg.e3=e3;

%ECEF spherical coordinates
xg.r=r; xg.theta=theta; xg.phi=phi;
xg.rx1i=[]; xg.thetax1i=[];
xg.rx2i=[]; xg.thetax2i=[];

%These are cartesian representations of the ECEF, spherical unit vectors
xg.er=er; xg.etheta=etheta; xg.ephi=ephi;

xg.I = p.I * ones([lx2,lx3]);

%Cartesian ECEF coordinates
xg.x=xECEF; xg.z=zECEF; xg.y=yECEF;
xg.alt=xg.r-Re;   %since we need a 3D array use xg.r here...

xg.gx1=gz; xg.gx2=zeros(xg.lx); xg.gx3=zeros(xg.lx);

xg.Bmag=-50000e-9*ones(xg.lx);     %minus for northern hemisphere...

%xg.glat = p.glat*ones(xg.lx);
xg.glon = p.glon*ones(xg.lx);    %use same lat./lon. for each grid point
xg.glat = glatgrid; xg.glon=glongrid;    %use same lat./lon. for each grid point

%xg.xp=x; xg.zp=z;

xg.inull=[];
%xg.nullpts=[];
xg.nullpts=zeros(xg.lx);


%% TRIM DATA STRUCTRE TO BE THE SIZE FORTRAN EXPECTS
xgf=xg;

% indices corresponding to non-ghost cells for 1 dimension
inds1=3:xgf.lx(1)-2;
inds2=3:xgf.lx(2)-2;
inds3=3:xgf.lx(3)-2;

% any dx variable will not need to first element (backward diff of two ghost cells)
indsdx1=2:xgf.lx(1);
indsdx2=2:xgf.lx(2);
indsdx3=2:xgf.lx(3);

% x1-interface variables need only non-ghost cell values (left interface) plus one
indsx1i=3:xgf.lx(1)-1;
indsx2i=3:xgf.lx(2)-1;
indsx3i=3:xgf.lx(3)-1;

% remove ghost cells
% now that indices have been define we can go ahead and make this change
xgf.lx=xgf.lx-4;

xgf.dx1b=xgf.dx1b(indsdx1);
xgf.dx2b=xgf.dx2b(indsdx2);
xgf.dx3b=xgf.dx3b(indsdx3);

xgf.x1i=xgf.x1i(indsx1i);
xgf.x2i=xgf.x2i(indsx2i);
xgf.x3i=xgf.x3i(indsx3i);

xgf.dx1h=xgf.dx1h(inds1);
xgf.dx2h=xgf.dx2h(inds2);
xgf.dx3h=xgf.dx3h(inds3);

xgf.h1x1i=xgf.h1x1i(indsx1i,inds2,inds3);
xgf.h2x1i=xgf.h2x1i(indsx1i,inds2,inds3);
xgf.h3x1i=xgf.h3x1i(indsx1i,inds2,inds3);

xgf.h1x2i=xgf.h1x2i(inds1,indsx2i,inds3);
xgf.h2x2i=xgf.h2x2i(inds1,indsx2i,inds3);
xgf.h3x2i=xgf.h3x2i(inds1,indsx2i,inds3);

xgf.h1x3i=xgf.h1x3i(inds1,inds2,indsx3i);
xgf.h2x3i=xgf.h2x3i(inds1,inds2,indsx3i);
xgf.h3x3i=xgf.h3x3i(inds1,inds2,indsx3i);

xgf.gx1=xgf.gx1(inds1,inds2,inds3);
xgf.gx2=xgf.gx2(inds1,inds2,inds3);
xgf.gx3=xgf.gx3(inds1,inds2,inds3);

xgf.glat=xgf.glat(inds1,inds2,inds3);
xgf.glon=xgf.glon(inds1,inds2,inds3);
xgf.alt=xgf.alt(inds1,inds2,inds3);

xgf.Bmag=xgf.Bmag(inds1,inds2,inds3);

xgf.I=xgf.I(inds2,inds3);

xgf.nullpts=xgf.nullpts(inds1,inds2,inds3);

xgf.e1=xgf.e1(inds1,inds2,inds3,:);
xgf.e2=xgf.e2(inds1,inds2,inds3,:);
xgf.e3=xgf.e3(inds1,inds2,inds3,:);

xgf.er=xgf.er(inds1,inds2,inds3,:);
xgf.etheta=xgf.etheta(inds1,inds2,inds3,:);
xgf.ephi=xgf.ephi(inds1,inds2,inds3,:);

xgf.r=xgf.r(inds1,inds2,inds3);
xgf.theta=xgf.theta(inds1,inds2,inds3);
xgf.phi=xgf.phi(inds1,inds2,inds3);

xgf.x=xgf.x(inds1,inds2,inds3);
xgf.y=xgf.y(inds1,inds2,inds3);
xgf.z=xgf.z(inds1,inds2,inds3);

end
