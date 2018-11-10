%LOWRES 2D EXAMPLE FOR TESTING
xdist=200e3;    %eastward distance
ydist=600e3;    %northward distance
lxp=80;
lyp=1;
glat=67.11;
glon=212.95;
gridflag=0;
I=90;


%ADD PATHS FOR FUNCTIONS
cwd = fileparts(mfilename('fullpath'));
addpath([cwd, filesep, '..', filesep,'..',filesep,'script_utils']);
addpath([cwd, filesep, '..', filesep,'..',filesep,'setup']);
addpath([cwd, filesep, '..', filesep,'..',filesep,'setup',filesep,'gridgen'])
addpath([cwd, filesep, '..', filesep,'..',filesep,'vis']);


%RUN THE GRID GENERATION CODE
if (~exist('xg'))
  xg=makegrid_cart_3D(xdist,lxp,ydist,lyp,I,glat,glon);
end
lx1=xg.lx(1); lx2=xg.lx(2); lx3=xg.lx(3);


%IDENTIFICATION FOR THE NEW SIMULATION THAT IS TO BE DONE
simid='2Dtest'


%ALTERNATIVELY WE MAY WANT TO READ IN AN EXISTING OUTPUT FILE AND DO SOME INTERPOLATION ONTO A NEW GRID
fprintf('Reading in source file...\n');
ID='../../../simulations/2Dtest_eq/'


%READ IN THE SIMULATION INFORMATION
[ymd0,UTsec0,tdur,dtout,flagoutput,mloc]=readconfig([ID,'/inputs/config.ini']);
xgin=readgrid([ID,'/inputs/']);


%FIND THE DATE OF THE END FRAME OF THE SIMULATION (PRESUMABLY THIS WILL BE THE STARTING POITN FOR ANOTEHR)
[ymdend,UTsecend]=dateinc(tdur,ymd0,UTsec0);


%LOAD THE FRAME
direc=ID;
[ne,mlatsrc,mlonsrc,xg,v1,Ti,Te,J1,v2,v3,J2,J3,filename,Phitop,ns,vs1,Ts]= ...
    loadframe(direc,ymdend,UTsecend,ymd0,UTsec0,tdur,dtout,flagoutput,mloc,xg);


%DO THE INTERPOLATION
[nsi,vs1i,Tsi]=model_resample(xgin,ns,vs1,Ts,xg);


%WRITE OUT THE GRID
outdir='../../../simulations/input/2Dtest/';
if (~(exist(outdir,'dir')==7))
  mkdir(outdir);
end
writegrid(xg,outdir);
dmy=[ymdend(3),ymdend(2),ymdend(1)];
writedata(dmy,UTsecend,nsi,vs1i,Tsi,outdir,simid);