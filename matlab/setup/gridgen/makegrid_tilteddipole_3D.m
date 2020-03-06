function xgf = makegrid_tilteddipole_3D(dtheta,dphi,lpp,lqp,lphip,altmin,glat,glon,gridflag)

narginchk(9,9)
%NOTE THAT INPUTS DTHETA AND DPHI ARE INTENDED TO REPRESENT THE FULL THETA
%AND PHI EXTENTS OF

%KNOWN ISSUES
%(1) For low latitude simulations your max L-shell should be large enough
%       that the altmin value is findable.  If it is not working increase pmax
%       to ~1.5 or so.
%(2) There is currently no way to produce a mesh which encompasses the
%       magnetic north pole
%(3) In plotting the cell structure the interface lines aren't truncated
%   properly.
%(5) Terrestrial mag. moment and mass are hard-coded in.
%(6) We assume target latitude is northern hemisphere
%
%OTHER NOTES:
% - Fortran code will interpret this grid as including ghost cells, so
% if you want a dimension to be size "n" adjust requested grid size so that
% it is "n+4"


%% PAD GRID WITH GHOST CELLS
lq=lqp+4;
lp=lpp+4;
lphi=lphip+4;


%% DEFINE DIPOLE GRID IN Q,P COORDS.
fprintf('\nMAKEGRID_TILTEDDIPOLE_3D.M --> Setting up q,p,phi grid of size %d x %d x %d.',lq-4,lp-4,lphi-4);
Re=6370e3;


%TD SPHERICAL LOCATION OF REQUESTED CENTER POINT
[thetatd,phid]=geog2geomag(glat,glon);


%SETS THE EDGES OF THE GRID
thetax2min=thetatd-dtheta/2*pi/180;
thetax2max=thetatd+dtheta/2*pi/180;
if(thetatd<pi/2)   %NH
  pmax=(Re+altmin)/Re/sin(thetax2min)^2;	%bottom left grid point p
  qtmp=(Re/(Re+altmin))^2*cos(thetax2min);	%bottom left grid q (also bottom right)
  pmin=sqrt(cos(thetax2max)/sin(thetax2max)^4/qtmp); %bottom right grid p
else               %SH
  pmax=(Re+altmin)/Re/sin(thetax2max)^2;	%bottom left grid point p
  qtmp=(Re/(Re+altmin))^2*cos(thetax2max);	%bottom left grid q (also bottom right)
  pmin=sqrt(cos(thetax2max)/sin(thetax2min)^4/qtmp); %bottom right grid p, why mixing of max/min here???
end
rtmp=fminbnd(@(x) qp2robj(x,qtmp,pmin),0,100*Re);        %bottom right r



%pmin=(Re+rtmp)/Re/sin(thetax2max)^2;
%p=linspace(pmin,pmax,lp);
p=linspace(pmin,pmax,lpp);
%p=p(:)';    %ensure a row vector
%pstride=p(2)-p(1);
%p=[p(1)-2*pstride,p(1)-pstride,p,p(end)+pstride,p(end)+2*pstride];

if gridflag==0      %open dipole grid
%    thetamax=thetax2min+pi/180;        %open
%    thetamax=thetax2min+pi/75;        %open
%     thetamax=thetamin+pi/50;        %open
%     thetamax=thetamin+pi/30;        %open
   if(thetatd<pi/2)   %northern hemisphere
     thetamax=thetax2min+pi/25;
   else
     thetamax=thetax2max-pi/25;
   end
else                %close dipole grid
   if(thetatd<pi/2) %NH
     thetamax=pi-thetax2min;
   else             %SH
     thetamax=pi-thetax2max;
   end
end
if(thetatd<pi/2)
  rmin=p(end)*Re*sin(thetax2min)^2; %use last field line to get qmin and qmax
  rmax=p(end)*Re*sin(thetamax)^2;
  qmin=cos(thetax2min)*Re^2/rmin^2;
  qmax=cos(thetamax)*Re^2/rmax^2;
else
  rmin=p(end)*Re*sin(thetamax)^2; %use last field line to get qmin and qmax
  rmax=p(end)*Re*sin(thetax2max)^2;
  qmin=cos(thetamax)*Re^2/rmin^2;
  qmax=cos(thetax2max)*Re^2/rmax^2;
end


%q=linspace(qmin,qmax,lq)';
q=linspace(qmin,qmax,lqp)';
q=sort(q);

p=p(:)';    %ensure a row vector
pstride=p(2)-p(1);
p=[p(1)-2*pstride,p(1)-pstride,p,p(end)+pstride,p(end)+2*pstride];

q=q(:);    %ensure a colume vector
qstride=q(2)-q(1);
q=[q(1)-2*qstride;q(1)-qstride;q,;q(end)+qstride;q(end)+2*qstride];    %add in ghost cells


%NOW THE AZIMUTHAL COORDINATE
phimin=phid-dphi/2*pi/180;
phimax=phid+dphi/2*pi/180;
%phi=linspace(phimin,phimax,lphi);    %note conversion to radians in  dphi calculation above
phi=linspace(phimin,phimax,lphip);
phi=phi(:)';
if (lphip>1)
  phistride=phi(2)-phi(1);     %assume constant stride
else
  phistride=0.1;   %just make up some junk for a 2D sim
end
phi=[phi(1)-2*phistride,phi(1)-phistride,phi,phi(end)+phistride,phi(end)+2*phistride];    %add in the ghost cells


%ALLOC/INIT - NOTE THESE DO NOT YET HAVE A PHI DIMENSION...
r=zeros(lq,lp);
fval=zeros(lp+1,lq+1);
theta=zeros(lq,lp);
qtol=1e-9;


%SPHERICAL XFORMATION
disp('MAKEGRID_TILTEDDIPOLE_3D: Converting q,p grid centers to spherical coords.')
for iq=1:lq
    for ip=1:lp
        [r(iq,ip),fval(iq,ip)]=fminbnd(@(x) qp2robj(x,q(iq),p(ip)),0,100*Re);
        if abs(q(iq))<qtol
            theta(iq,ip)=pi/2;
        elseif q(iq)>=0        %northern hemisphere
            theta(iq,ip)=asin(sqrt(r(iq,ip)/Re/p(ip)));
        else
            theta(iq,ip)=pi-asin(sqrt(r(iq,ip)/Re/p(ip)));
        end
    end
end


%EXTEND R,TH ARRAYS TO BE 3D
r=repmat(r,[1 1 lphi]);
theta=repmat(theta,[1 1 lphi]);
phispher=repmat(reshape(phi(:),[1 1 lphi]),[lq,lp,1]);


%TRUE CARTESIAN
z=r.*cos(theta);
x=r.*sin(theta).*cos(phispher);
y=r.*sin(theta).*sin(phispher);


%CARTESIAN FOR PLOTTING PURPOSES
%{
r1=mean(r(1,:));
r2=mean(r(lq,:));
if gridflag==0
    if r1<r2
        meanth=mean(theta(1,:));
    else
        meanth=mean(theta(lq,:));
    end
else
    meanth=pi/2;
end
xsurf=Re*sin(meanth);
zsurf=Re*cos(meanth);
xp=x-xsurf;
zp=z-zsurf;
rotmat=[cos(meanth), -sin(meanth); sin(meanth), cos(meanth)];
for iq=1:lq
    for ip=1:lp
        xz=rotmat*[xp(iq,ip); zp(iq,ip)];
        xp(iq,ip)=xz(1); zp(iq,ip)=xz(2);
    end
end
% x=r.*sin(theta-meanth);
% z=r.*cos(theta-meanth);
%}


%INTERFACE LOCATIONS
disp('MAKEGRID_TILTEDDIPOLE_3D.M --> Converting q,p grid interfaces to spherical coords.');
qi=zeros(lq+1,1);
qi(2:lq)=1/2*(q(1:lq-1)+q(2:lq));
qi(1)=q(1)-1/2*(q(2)-q(1));
qi(lq+1)=q(lq)+1/2*(q(lq)-q(lq-1));
pii=zeros(1,lp+1);
pii(2:lp)=1/2*(p(1:lp-1)+p(2:lp));
pii(1)=p(1)-1/2*(p(2)-p(1));
pii(lp+1)=p(lp)+1/2*(p(lp)-p(lp-1));
phii=zeros(lphi+1,1);
phii(2:lphi)=1/2*(phi(1:lphi-1)+phi(2:lphi));
phii(1)=phi(1)-1/2*(phi(2)-phi(1));
phii(lphi+1)=phi(lphi)+1/2*(phi(lphi)-phi(lphi-1));


%SPHERICAL XFORMATION OF INTERFACES IN Q
rqi=zeros(lq+1,lp);
thetaqi=zeros(lq+1,lp);
for iq=1:lq+1
    for ip=1:lp
        [rqi(iq,ip),fval(iq,ip)]=fminbnd(@(x) qp2robj(x,qi(iq),p(ip)),0,100*Re);
        if abs(qi(iq))<qtol
            thetaqi(iq,ip)=pi/2;
        elseif qi(iq)>=0        %northern hemisphere
            thetaqi(iq,ip)=asin(sqrt(rqi(iq,ip)/Re/p(ip)));
        else
            thetaqi(iq,ip)=pi-asin(sqrt(rqi(iq,ip)/Re/p(ip)));
        end
    end
end


%COPY QI INFO FOR 3RD DIM
rqi=repmat(rqi,[1 1 lphi]);
thetaqi=repmat(thetaqi,[1 1 lphi]);


%SPHERICAL XFORMATION OF INTERFACES IN P
rpi=zeros(lq,lp+1);
thetapi=zeros(lq,lp+1);
for iq=1:lq
    for ip=1:lp+1
        [rpi(iq,ip),fval(iq,ip)]=fminbnd(@(x) qp2robj(x,q(iq),pii(ip)),0,100*Re);
        if abs(q(iq))<qtol
            thetapi(iq,ip)=pi/2;
        elseif q(iq)>=0        %northern hemisphere
            thetapi(iq,ip)=asin(sqrt(rpi(iq,ip)/Re/pii(ip)));
        else
            thetapi(iq,ip)=pi-asin(sqrt(rpi(iq,ip)/Re/pii(ip)));
        end
    end
end


%COPY PI INFO FOR 3RD DIM
rpi=repmat(rpi,[1 1 lphi]);
thetapi=repmat(thetapi,[1 1 lphi]);


%METRIC COEFFICIENTS
disp('MAKEGRID_TILTEDDIPOLE_3D.M --> Calculating metric coeffs.');
denom=sqrt(1+3*cos(theta).^2);
hq=r.^3/Re^2./denom;
hp=Re*sin(theta).^3./denom;
hphi=r.*sin(theta);

hqphii=cat(3,hq,hq(:,:,lphi)); %phi interfaces at exactly r,th coords of cell centers.
hpphii=cat(3,hp,hp(:,:,lphi));
hphiphii=cat(3,hphi,hphi(:,:,lphi));

denom=sqrt(1+3*cos(thetaqi).^2);
hqqi=rqi.^3/Re^2./denom;
hpqi=Re*sin(thetaqi).^3./denom;
hphiqi=rqi.*sin(thetaqi);

denom=sqrt(1+3*cos(thetapi).^2);
hqpi=rpi.^3/Re^2./denom;
hppi=Re*sin(thetapi).^3./denom;
hphipi=rpi.*sin(thetapi);


%SPHERICAL UNIT VECTORS IN CARTESIAN COMPONENTS (CELL-CENTERED)
disp('MAKEGRID_TILTEDDIPOLE_3D.M --> Calculating spherical unit vectors.');
er(:,:,:,1)=sin(theta).*cos(phispher);
er(:,:,:,2)=sin(theta).*sin(phispher);
er(:,:,:,3)=cos(theta);
etheta(:,:,:,1)=cos(theta).*cos(phispher);
etheta(:,:,:,2)=cos(theta).*sin(phispher);
etheta(:,:,:,3)=-sin(theta);
ephi(:,:,:,1)=-sin(phispher);
ephi(:,:,:,2)=cos(phispher);
ephi(:,:,:,3)=zeros(lq,lp,lphi);


%UNIT VECTORS FOR Q,P,PHI FOR ALL GRID POINTS IN CARTESIAN COMPONENTS
disp('MAKEGRID_TILTEDDIPOLE_3D.M --> Calculating dipole unit vectors.');
denom=Re^2*(1+3*cos(theta).^2);
dxdq(:,:,:,1)=-3*r.^3.*cos(theta).*sin(theta)./denom.*cos(phispher);
dxdq(:,:,:,2)=-3*r.^3.*cos(theta).*sin(theta)./denom.*sin(phispher);
dxdq(:,:,:,3)=(-2*cos(theta).^2+sin(theta).^2).*r.^3./denom;
magdxdq=repmat(sqrt(dot(dxdq,dxdq,4)),[1,1,1,3]);
eq=dxdq./magdxdq;
ep=cross(ephi,eq,4);
Imat=acos(dot(er,eq,4));
if gridflag==0
    I=mean(Imat,1);             %avg. inclination for each field line.
else
    I=mean(Imat(1:floor(lq/2),:,:),1);   %avg. over only half the field line
end
I=90-min(I,pi-I)*180/pi;    %ignore parallel vs. anti-parallel


%GRAVITATIONAL FIELD COMPONENTS IN DIPOLE SYSTEM
G=6.67428e-11;
Me=5.9722e24;
g=G*Me./r.^2;
gq=g.*dot(-er,eq,4);
gp=g.*dot(-er,ep,4);


%MAGNETIC FIELD STRENGTH
Bmag=(4*pi*1e-7)*7.94e22/4/pi./(r.^3).*sqrt(3*(cos(theta)).^2+1);


%STORE RESULTS IN GRID DATA STRUCTURE
disp('MAKEGRID_TILTEDDIPOLE_3D.M --> Creating a grid structure with the results.\n');
xg.x1=q; xg.x2=p; xg.x3=reshape(phi,[1 1 lphi]);
xg.x1i=qi; xg.x2i=pii; xg.x3i=reshape(phii,[1 1 lphi+1]);
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

xg.dx3f=zeros(1,1,lphi);
xg.dx3b=xg.dx3f;
xg.dx3h=xg.dx3f;
dx1=xg.x3(2)-xg.x3(1);
dxn=xg.x3(lx(3))-xg.x3(lx(3)-1);
xg.dx3f=cat(3,xg.x3(1,1,2:lx(3))-xg.x3(1,1,1:lx(3)-1), dxn);         %FWD DIFF
xg.dx3b=cat(3,dx1, xg.x3(1,1,2:lx(3))-xg.x3(1,1,1:lx(3)-1));         %BACK DIFF
xg.dx3h=xg.x3i(1,1,2:lx(3)+1)-xg.x3i(1,1,1:lx(3));              %MIDPOINT DIFFS

xg.h1=hq; xg.h2=hp; xg.h3=hphi;
xg.h1x1i=hqqi; xg.h2x1i=hpqi; xg.h3x1i=hphiqi;
xg.h1x2i=hqpi; xg.h2x2i=hppi; xg.h3x2i=hphipi;
xg.h1x3i=hqphii; xg.h2x3i=hpphii; xg.h3x3i=hphiphii;

xg.e1=eq; xg.e2=ep; xg.e3=ephi;

xg.r=r; xg.theta=theta; xg.phi=phispher;
xg.rx1i=rqi; xg.thetax1i=thetaqi;
xg.rx2i=rpi; xg.thetax2i=thetapi;

xg.er=er; xg.etheta=etheta; xg.ephi=ephi;

xg.I=I;

xg.x=x; xg.z=z; xg.y=y;
xg.alt=xg.r-Re;

%xg.xp=xp; xg.zp=zp;

xg.gx1=gq; xg.gx2=gp; xg.gx3=zeros(lq,lp,lphi);

xg.Bmag=Bmag;

%xg.glat=(pi/2-theta)*180/pi; xg.glon=phi*180/pi*ones(lx(1),lx(2));
for iphi=1:lphi
  [glats,glons]=geomag2geog(xg.theta(:,:,iphi),xg.phi(1,1,iphi)*ones(lq,lp));    %only meant to work for one latitude at at time
  xg.glat(:,:,iphi)=glats;
  xg.glon(:,:,iphi)=glons;
end

% xg.inull=find(r<Re+30e3); %may give issues in conservative form???  NOPE not the problem
% xg.nullpts=r<Re+30e3;
xg.inull=find(r<Re+80e3);
xg.nullpts=r<Re+80e3;

%NOW ADJUST SIZES SO THAT THEY MATCH WHAT FORTRAN CODE EXPECTS.  IF NOT
%USING THIS TO GENERATE A GRID FOR THE FORTRAN CODE YOU MAY WANT TO GET RID
%OF THIS BLOCK
xgf=xg;    %make a copy to alter for purposes of getting all of the sizes the same as used by fortran

inds1=3:xgf.lx(1)-2;    %indices corresponding to non-ghost cells for 1 dimension
inds2=3:xgf.lx(2)-2;
inds3=3:xgf.lx(3)-2;

indsdx1=2:xgf.lx(1);    %any dx variable will not need to first element (backward diff of two ghost cells)
indsdx2=2:xgf.lx(2);
indsdx3=2:xgf.lx(3);

indsx1i=3:xgf.lx(1)-1;    %x1-interface variables need only non-ghost cell values (left interface) plus one
indsx2i=3:xgf.lx(2)-1;
indsx3i=3:xgf.lx(3)-1;

xgf.lx=xgf.lx-4;    %remove ghost cells, now that indices have been define we can go ahead and make this change

xgf.dx1b=xgf.dx1b(indsdx1);
xgf.dx2b=xgf.dx2b(indsdx2);
xgf.dx3b=xgf.dx3b(indsdx3);

xgf.x1i=xgf.x1i(indsx1i);
xgf.x2i=xgf.x2i(indsx2i);
xgf.x3i=xgf.x3i(indsx3i);

xgf.dx1h=xgf.dx1h(inds1);
xgf.dx2h=xgf.dx2h(inds2);
xgf.dx3h=xgf.dx3h(inds3);

xgf.h1=xgf.h1;    %these need ghost cells for compression divergence term
xgf.h2=xgf.h2;
xgf.h3=xgf.h3;

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

xgf.I=xgf.I(1,inds2,inds3);

xgf.nullpts=xgf.nullpts(inds1,inds2,inds3);

%ZZZ - NEED TO ALSO CORRECT OTHER VARIABLE SIZES!!!!
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
