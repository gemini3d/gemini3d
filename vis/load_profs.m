function [xg,neprof,Teprof,Tiprof,viprof,simdate] = load_profs(direc,xg)


%SET PATHS FOR FUNCTIONS NEEDED INSIDE THIS ROUTINE
cwd = fileparts(mfilename('fullpath'));
addpath([cwd, filesep, 'plotfunctions'])
addpath([cwd, filesep, '..', filesep, 'script_utils'])


%%READ IN THE SIMULATION INFORMATION
[ymd0,UTsec0,tdur,dtout,flagoutput,mloc]=readconfig([direc,filesep,'inputs/config.ini']);


%CHECK WHETHER WE NEED TO RELOAD THE GRID (WHICH CAN BE TIME CONSUMING)
if (~exist('xg','var'))
  xg = readgrid([direc,filesep,'inputs',filesep]);
end


% COMPUTE SOURUCE LOCATION IN MCOORDS
if (~isempty(mloc))
 mlat=mloc(1);
 mlon=mloc(2);
else
 mlat=[];
 mlon=[];
end


%% TIMES OF INTEREST
times=UTsec0:dtout:UTsec0+tdur;
lt=numel(times);


%LOCATIONS OF THE FIELD ALIGNED PROFILES TO BE EXTRACTED
ix2=min(find(xg.x2(3:end-2)>0));
ix3=min(find(xg.x3(3:end-2)>0))+floor(xg.lx(3)/8);


%ALLOCATE SPACE FOR TIME SERIES
neprof=zeros(xg.lx(1),lt);
Tiprof=zeros(xg.lx(1),lt);
Teprof=zeros(xg.lx(1),lt);
viprof=zeros(xg.lx(1),lt);


%MAIN TIME LOOP FOR LOADING AND EXTRACTING DATA OF INTEREST
ymd=ymd0;
UTsec=UTsec0;
simdate=[];
for it=1:lt
    simdate=cat(1,simdate,[ymd,UTsec/3600,0,0]);
    [ne,mlatsrc,mlonsrc,v1,Ti,Te,J1,v2,v3,J2,J3,filename,Phitop] = loadframe(direc,UTsec,ymd,UTsec0,ymd0,mloc,xg);
    disp(filename)

    neprof(:,it)=ne(:,ix2,ix3);
    Teprof(:,it)=Te(:,ix2,ix3);
    Tiprof(:,it)=Ti(:,ix2,ix3);
    viprof(:,it)=v1(:,ix2,ix3);
    
    [ymd,UTsec]=dateinc(dtout,ymd,UTsec);
end % for
    
    
end % function
