%SIMULATIONS LOCAITONS
simname='test3d_fang_mag/';
basedir='../../../simulations/'
direc=[basedir,simname];


%SIMULATION META-DATA
[ymd0,UTsec0,tdur,dtout,flagoutput,mloc]=readconfig([direc,'/inputs']);


%ANGULAR RANGE TO COVER FOR THE CALCLUATIONS (THIS IS FOR THE FIELD POINTS - SOURCE POINTS COVER ENTIRE GRID)
dang=1.5;


%WE ALSO NEED TO LOAD THE GRID FILE
if ~exist('xg','var')
  disp('Reading grid...')
  xg=readgrid([direc,'/inputs/']);
  lx1=xg.lx(1); lx2=xg.lx(2); lx3=xg.lx(3);
  lh=lx1;   %possibly obviated in this version - need to check
  if (lx3==1)
    flag2D=1;
    fprintf('2D meshgrid...\n')
    x1=xg.x1(3:end-2);
    x2=xg.x2(3:end-2);
    x3=xg.x3(3:end-2);
    [X2,X1]=meshgrid(x2(:),x1(1:lh)');
  else
    flag2D=0;
    fprintf('3D meshgrid...\n')
    x1=xg.x1(3:end-2);
    x2=xg.x2(3:end-2);
    x3=xg.x3(3:end-2);
    [X2,X1,X3]=meshgrid(x2(:),x1(1:lh)',x3(:));
  end


  %TABULATE THE SOURCE OR GRID CENTER LOCATION
  if (~isempty(mloc))
    mlatsrc=mloc(1);
    mlonsrc=mloc(2);
    thdist=pi/2-mlatsrc*pi/180;    %zenith angle of source location
    phidist=mlonsrc*pi/180;
  else
    thdist=mean(xg.theta(:));
    phidist=mean(xg.phi(:));
  end
end
fprintf('Grid loaded...\n');


%FIELD POINTS OF INTEREST (CAN/SHOULD BE DEFINED INDEPENDENT OF SIMULATION GRID)
ltheta=10;
if (~flag2D)
  lphi=10;
else
  lphi=1;
end
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
r=6370e3*ones(ltheta,lphi);                          %use ground level for altitude for all field points
[phi,theta]=meshgrid(phi,theta);

%CREATE AN INPUT FILE OF FIELD POINTS
fid=fopen([basedir,'/input/',simname,'/magfieldpoints.dat'],'w');
fwrite(fid,numel(theta),'integer*4');
fwrite(fid,r(:),'real*8');
fwrite(fid,theta(:),'real*8');
fwrite(fid,phi(:),'real*8');
