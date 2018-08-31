function xg=makegrid_tilteddipole(dtheta,dphi,lp,lq,lphi,altmin,glat,glon,gridflag,plotflag)

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


%DEFINE DIPOLE GRID IN Q,P COORDS.
fprintf('\nMAKEGRID_DIPOLE.M --> Setting up q,p grid of size %d x %d.',lq,lp);
Re=6370e3;


%TD SPHERICAL LOCATION OF REQUESTED CENTER POINT
[thetatd,phid]=geog2geomag(glat,glon);

thetax2min=thetatd-dtheta/2*pi/180;
thetax2max=thetatd+dtheta/2*pi/180;
pmax=(Re+altmin)/Re/sin(thetax2min)^2;	%bottom left grid point p
qtmp=(Re/(Re+altmin))^2*cos(thetax2min);	%bottom left grid q (also bottom right)
pmin=sqrt(cos(thetax2max)/sin(thetax2max)^4/qtmp); %bottom right grid p
rtmp=fminbnd(@(x) qp2robj(x,qtmp,pmin),0,100*Re);        %bottom right r
%pmin=(Re+rtmp)/Re/sin(thetax2max)^2;
p=linspace(pmin,pmax,lp);
if gridflag==0
%    thetamax=thetamin+pi/100;        %open
%    thetamax=thetamin+pi/75;        %open
%     thetamax=thetamin+pi/50;        %open
%     thetamax=thetamin+pi/30;        %open
    thetamax=thetax2min+pi/25;        %open
else
    thetamax=pi-thetax2min;           %closed
end
rmin=p(lp)*Re*sin(thetax2min)^2; %use last field line to get qmin and qmax
rmax=p(lp)*Re*sin(thetamax)^2;
qmin=cos(thetax2min)*Re^2/rmin^2;
qmax=cos(thetamax)*Re^2/rmax^2;
q=linspace(qmin,qmax,lq)';
q=sort(q);


%NOW THE AZIMUTHAL COORDINATE
phimin=phid-dphi/2;
phimax=phid+dphi/2;
phi=linspace(phimin*pi/180,phimax*pi/180,lphi);    %note conversion to radians


%ALLOC/INIT - NOTE THESE DO NOT YET HAVE A PHI DIMENSION...
r=zeros(lq,lp);
fval=zeros(lp+1,lq+1);
theta=zeros(lq,lp);
qtol=1e-9;


%SPHERICAL XFORMATION
fprintf('\nMAKEGRID_DIPOLE.M --> Converting q,p grid centers to spherical coords.');
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
y=r.*sin(theta).*sin(phispher)


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
fprintf('\nMAKEGRID_DIPOLE.M --> Converting q,p grid interfaces to spherical coords.');
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
fprintf('\nMAKEGRID_DIPOLE.M --> Calculating metric coeffs. and unit vectors.\n');
denom=sqrt(1+3*cos(theta).^2);
hq=r.^3/Re^2./denom;
hp=Re*sin(theta).^3./denom;
hphi=r.*sin(theta);

denom=sqrt(1+3*cos(thetaqi).^2);
hqqi=rqi.^3/Re^2./denom;
hpqi=Re*sin(thetaqi).^3./denom;
hphiqi=rqi.*sin(thetaqi);

denom=sqrt(1+3*cos(thetapi).^2);
hqpi=rpi.^3/Re^2./denom;
hppi=Re*sin(thetapi).^3./denom;
hphipi=rpi.*sin(thetapi);


%SPHERICAL UNIT VECTORS IN CARTESIAN COMPONENTS (CELL-CENTERED)
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
denom=Re^2*(1+3*cos(theta).^2);
dxdq(:,:,:,1)=-3*r.^3.*cos(theta).*sin(theta)./denom.*cos(phispher);
dxdq(:,:,:,2)=-3*r.^3.*cos(theta).*sin(theta)./denom.*sin(phispher);
dxdq(:,:,:,3)=(-2*cos(theta).^2+sin(theta).^2).*r.^3./denom;
magdxdq=repmat(sqrt(dot(dxdq,dxdq,4)),[1,1,1,3]);
eq=dxdq./magdxdq;
ep=cross(ephi,eq,4);
Imat=acos(dot(er,eq,4));
if gridflag==0    I=mean(Imat,1);             %avg. inclination for each field line.
else
    I=mean(Imat(1:floor(lq/2),:,:),1);   %avg. over only half the field line
end
I=90-min(I,pi-I)*180/pi;    %ignore parallel vs. anti-parallel


%DIAGNOSTIC PLOTS
if plotflag
    figure;
    set(gcf,'PaperPosition',[0 0 8.5 4]);
    subplot(121);
    polar(theta(:),r(:),'k.');
    hold on;
    %     polar(thetaqi(:),rqi(:),'r.');
    %     polar(thetapi(:),rpi(:),'g.');
    for iq=1:lq+1
        polar(thetaqi(iq,:),rqi(iq,:));
    end
    for ip=1:lp+1
        polar(thetapi(:,ip),rpi(:,ip));
    end
    hold off;
    
    subplot(122);
    plot(x,z,'.')
    hold on;
    % quiver(x,z,er(:,:,1),er(:,:,3));
    % quiver(x,z,etheta(:,:,1),etheta(:,:,3));
    %     quiver(x,z,eq(:,:,1),eq(:,:,3));
    %     quiver(x,z,ep(:,:,1),ep(:,:,3));
    hold off;
end


%GRAVITATIONAL FIELD COMPONENTS IN DIPOLE SYSTEM
G=6.67428e-11;
Me=5.9722e24;
g=G*Me./r.^2;
gq=g.*dot(-er,eq,4);
gp=g.*dot(-er,ep,4);


%MAGNETIC FIELD STRENGTH
Bmag=(4*pi*1e-7)*7.94e22/4/pi./(r.^3).*sqrt(3*(cos(theta)).^2+1);


%STORE RESULTS IN GRID DATA STRUCTURE
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

xg.r=r; xg.theta=theta; xg.phi=phi;
xg.rx1i=rqi; xg.thetax1i=thetaqi;
xg.rx2i=rpi; xg.thetax2i=thetapi;

xg.er=er; xg.etheta=etheta; xg.ephi=ephi;

xg.I=I;

xg.x=x; xg.z=z; xg.y=y;
xg.alt=xg.r-Re;

xg.xp=xp; xg.zp=zp;

xg.gx1=gq; xg.gx2=gp; xg.gx3=zeros(lq,lp,lphi);

xg.Bmag=Bmag;

%xg.glat=(pi/2-theta)*180/pi; xg.glon=phi*180/pi*ones(lx(1),lx(2));
[glats,glons]=geomag2geog(xg.theta,xg.phi);
xg.glat=glats;
xg.glon=glons;

% xg.inull=find(r<Re+30e3); %may give issues in conservative form???  NOPE not the problem
% xg.nullpts=r<Re+30e3;
xg.inull=find(r<Re+80e3);
xg.nullpts=r<Re+80e3;

end
