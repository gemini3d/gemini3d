
%MAKE A GRID
% dtheta=25;
% dphi=45;
% lp=35;
% lq=350;
% lphi=35;
% altmin=80e3;
% glat=35;
% glon=270;
% gridflag=0;

dtheta=2;
dphi=5;
lp=110;
lq=200;
lphi=80;
altmin=80e3;
glat=65;    %high-latitude
glon=270;
gridflag=0;

xg=makegrid_tilteddipole_3D(dtheta,dphi,lp,lq,lphi,altmin,glat,glon,gridflag);


%GENERATE INITIAL CONDITIONS
fprintf('Generating intial conditions...\n');
UT=0;
dmy=[15,9,2016];
activ=[100,100,10];
nmf=5e11;
nme=2e11;
[ns,Ts,vsx1]=eqICs3D(xg,UT,dmy,activ,nmf,nme);    %note that this actually calls msis_matlab - should be rewritten to includ the neutral module form the fortran code!!!


%WRITE THE GRID AND INITIAL CONDITIONS
fprintf('Writing grid to file...\n');
writegrid(xg,'curvtest');

fprintf('Writing initial conditions to file...\n');
time=UT*3600;   %doesn't matter for input files
writedata(time,ns,Ts,vsx1,'curvtest');
